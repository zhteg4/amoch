// ============================================================================
// restart.cc -- Functions to restart a simulation
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#include <cstdio>
#include <cstring>
#include "restart.h"
#include "logfile.h"
#include "os.h"

using std::vector;
using AMOCH::System;
using AMOCH::Chain;
using AMOCH::log_message;

// ============================================================================
// Return the index within the Chain chains[jc] of particle id; assumes the
// LAMMPS id comes from the atoms.dat file in which Chains are written in 
// order
// ============================================================================
static int
convert_id_to_index(const vector<Chain>& chains,
                    int jc,
                    int id)
{
   int n = 0;

   for (int i = 0; i < jc; i++)
      n += chains[i].nparticles();
   return id - n;
}

// ============================================================================
// Return true if particles in LAMMPS dump file specified by path can be read
// into sys, else log an error message and return false.  This function expects
// a specific column order and UNwrapped coordinates; the dump file can be 
// written using these arguments with a custom dump: "id mol type xu yu zu"
// For example,
//    dump 1 all custom 1 file.dump id mol type xu yu zu
// ============================================================================
bool
AMOCH::read_lammps_dump(FILE *f,
                        System& sys,
                        vector<Chain>& chains)
{
   int N, i;
   const size_t linelen = 128;
   char line[linelen];

   // Number of atoms
   do {
      fgets(line, linelen, f);
   } while (!strstr(line, "ITEM: NUMBER OF ATOMS"));
   fgets(line, linelen, f);
   sscanf(line, "%d", &N);

   // ATOMS format
   while (fgets(line, linelen, f)) {
      if (strstr(line, "ITEM: ATOMS id mol type xu yu zu"))
         break;
   }
   if (feof(f) || ferror(f)) {
      log_message("\nDid not find expected \"ITEM: ATOMS\" format in ");
      log_message("LAMMPS dump file\n");
      return false;
   }

   // Particles
   int id, mol, type, index;
   int cell, slot;
   AMOCH::Triple p, tp;

   for (i = 0; i < N; i++) {
      if (!fgets(line, linelen, f))
         break;
      sscanf(line, "%d%d%d" SCN_REAL SCN_REAL SCN_REAL, &id, &mol, &type, &tp.x,
             &tp.y, &tp.z);
      --id;
      --mol;
      --type;
      index = convert_id_to_index(chains, mol, id);
      p = tp;
      sys.wrap(p);
      cell = sys.cell_index(p);
      slot = sys.cells[cell].add_particle(mol, index, type, p, tp);
      chains[mol].add_particle_cell(index, cell, slot);
   }
   if (i < N) {
      log_message("\nDid not read %d atoms from LAMMPS dump file\n", N);
      return false;
   }
   return true;
}

// ============================================================================
// Return the (Chain) index within chains of particle id
// ============================================================================
static int
find_chain(const vector<Chain>& chains,
           int nchains,
           int id)
{
   int n = 0;
   int i;

   for (i = 0; i < nchains-1; i++) {
      if (id < n+chains[i].nparticles())
         break;
      n += chains[i].nparticles();
   }
   return i;
}

// ============================================================================
// Return true if the bond info in the indicated file can be read into sys,
// else return false after logging a message
// ============================================================================
bool
AMOCH::read_bonds(FILE *f,
                  System& sys,
                  vector<Chain>& chains)
{
   bool retval = true;
   const size_t linelen = 128;
   char line[linelen];
   const int nchains = chains.size();
   int id1, id2, bid, btype;
   int mol, index1, index2, cell, slot;

   fgets(line, linelen, f);
   fgets(line, linelen, f);
   while (fgets(line, linelen, f)) {
      sscanf(line, "%d%d%d%d", &bid, &btype, &id1, &id2);
      --id1;
      --id2;
      --btype;
      mol = find_chain(chains, nchains, id1);  // same for id2
      index1 = convert_id_to_index(chains, mol, id1);
      index2 = convert_id_to_index(chains, mol, id2);
      cell = chains[mol].cell(index1);
      slot = chains[mol].slot(index1);
      sys.cells[cell].add_bond(slot, index2, btype);
      cell = chains[mol].cell(index2);
      slot = chains[mol].slot(index2);
      sys.cells[cell].add_bond(slot, index1, btype);
   }
   if (ferror(f)) {
      log_message("\nError reading bonds\n");
      retval = false;
   }
   return retval;
}

