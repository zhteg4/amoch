// ============================================================================
// output.cc -- Output functions
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include "output.h"

using std::vector;
using AMOCH::Real;
using AMOCH::System;
using AMOCH::Chain;
using AMOCH::PartType;

// ============================================================================
// Write a status bar with completion percentage to the output FILE f
// ============================================================================
void
AMOCH::show_status(FILE *f,
                   int step,
                   int nsteps)
{
   const int barlen = 68;
   const Real r = (Real)step / (Real)nsteps;
   const int n = (int)(barlen * r);

   fputc('|', f);
   for (int i = 0; i < barlen; i++)
      fputc(((i <= n) ? '=' : ' '), f);
   fprintf(f, "| %6.2f %%%c", 100.0*r, ((step == nsteps) ? '\n' : '\r'));
   fflush(f);
}
      
// ============================================================================
// Write a PDB output file of all Chains in a System
// ============================================================================
void
AMOCH::write_pdb(FILE *f, 
                 const char *comment,
                 const System& sys,
                 const vector<Chain>& chains,
                 const vector<PartType>& part_types,
                 bool wrapped,
                 bool connect)
{
   unsigned n = 1;
   unsigned i;
   AMOCH::Triple p;
   int j, k, cell, slot;
   const char *csym;
   const char *atom_fmt = 
      "ATOM  %5d %4s              %8.3f%8.3f%8.3f                      %2s  \n";

   fprintf(f, "REMARK   1\nREMARK   1 %s\n", comment);
   for (i = 0; i < chains.size(); i++) {
      for (j = 0; j < chains[i].nparticles(); j++, n++) {
         cell = chains[i].cell(j);
         slot = chains[i].slot(j);
         p = (wrapped) ? sys.cells[cell].pos[slot] : sys.cells[cell].tpos[slot];
         csym = part_types[sys.cells[cell].ptype(slot)].symbol().c_str();
         fprintf(f, atom_fmt, n, csym, p.x, p.y, p.z, csym);
      }
      fprintf(f, "TER   %5d\n", n++);
   }
   if (connect) {
      int offset = 0;

      n = 1;
      for (i = 0; i < chains.size(); i++) {
         for (j = 0; j < chains[i].nparticles(); j++, n++) {
            cell = chains[i].cell(j);
            slot = chains[i].slot(j);
            fprintf(f, "CONECT%5d", n);
            for (k = 0; k < sys.cells[cell].nbonds(slot); k++)
               fprintf(f, "%5d", offset+sys.cells[cell].bonds(slot)[k][0]+1);
            fprintf(f, "\n");
         }
         offset += j+1;
         n++;  // TER record -- no corresponding CONECT
      }
   }
   fprintf(f, "END\n");
   fflush(f);
}

// ============================================================================
// Write a XYZ output file of all Chains in a System
// ============================================================================
void
AMOCH::write_xyz(FILE *f, 
                 const char *comment,
                 const System& sys,
                 const vector<Chain>& chains,
                 const vector<PartType>& part_types,
                 bool wrapped)
{
   unsigned i;
   AMOCH::Triple p;
   int j, cell, slot;
   int n = 0;
   const char *csym;
   const char *fmt = "%4s%15.3f%15.3f%15.3f\n";

   for (i = 0; i < chains.size(); i++)
      n += chains[i].nparticles();
   fprintf(f, "%d\n%s\n", n, comment);
   for (i = 0; i < chains.size(); i++) {
      for (j = 0; j < chains[i].nparticles(); j++) {
         cell = chains[i].cell(j);
         slot = chains[i].slot(j);
         p = (wrapped) ? sys.cells[cell].pos[slot] : sys.cells[cell].tpos[slot];
         csym = part_types[sys.cells[cell].ptype(slot)].symbol().c_str();
         fprintf(f, fmt, csym, p.x, p.y, p.z);
      }
   }
   fflush(f);
}

// ============================================================================
// Write atoms.dat file for data4Lammps
// ============================================================================
void 
AMOCH::write_atoms_dat(FILE *f,
                       const System& sys, 
                       const vector<Chain>& chains)
{
   unsigned i;
   int n, j, cell, slot;

   n = 0;
   for (i = 0; i < chains.size(); i++)
      n += chains[i].nparticles();
   fprintf(f, "%d atoms\n\nATOMS\n\n", n);
   n = 1;
   for (i = 0; i < chains.size(); i++) {
      for (j = 0; j < chains[i].nparticles(); j++) {
         cell = chains[i].cell(j);
         slot = chains[i].slot(j);
         fprintf(f, "%d  %d  %d  %.6f  %.12f  %.12f  %.12f\n", 
                 n++, 
                 i+1, 
                 sys.cells[cell].ptype(slot)+1, 
                 0.0,  // XXX charge
                 sys.cells[cell].pos[slot].x,
                 sys.cells[cell].pos[slot].y,
                 sys.cells[cell].pos[slot].z);
      }
   }
   fflush(f);
}
                 
// ============================================================================
// Write atom_type.dat input file for data4Lammps; the particle types are not
// specific to any force field, but all the info is written for subsequent
// forcefield typing
// ============================================================================
void 
AMOCH::write_atom_type_dat_in(FILE *f, 
                              const vector<PartType>& part_types)
{
   unsigned i;
   int n;

   n = 1;
   for (i = 0; i < part_types.size(); i++) {
      fprintf(f, "%d  %s  %d  %d  %d\n", 
              n++,
              part_types[i].symbol().c_str(),  // element string
              part_types[i].nbonds(),
              part_types[i].resonant(),
              part_types[i].hbond());
   }
   fflush(f);
}

// ============================================================================
// Write bonds.dat input file for data4Lammps
// ============================================================================
void
AMOCH::write_bonds_dat(FILE *f, 
                       const System& sys, 
                       const vector<Chain>& chains)
{
   unsigned i;
   int j, k, n;
   int jndx, kndx;
   int offset = 0;
   int btype;
   int jcell, jslot;
   AMOCH::Ivec bond(2);

   fprintf(f, "BONDS\n\n");
   n = 1;
   for (i = 0; i < chains.size(); i++) {
      for (j = 0; j < chains[i].nparticles(); j++) {
         jcell = chains[i].cell(j);
         jslot = chains[i].slot(j);
         jndx = sys.cells[jcell].index(jslot);
         for (k = 0; k < sys.cells[jcell].nbonds(jslot); k++) {
            bond = sys.cells[jcell].bonds(jslot)[k];
            kndx = bond[0];
            btype = bond[1];
            if (jndx < kndx)  {  // write bond once
               fprintf(f, "%d  %d  %d  %d\n", 
                       n++, 
                       btype+1,
                       offset+jndx+1,
                       offset+kndx+1);
            }
         }
      }
      offset += j;
   }
   fflush(f);
}

// ============================================================================
// Write bond_type.dat input file for data4Lammps; the atom types are written
// as indices for types in atom_type.dat, not specific to any force field
// ============================================================================
void 
AMOCH::write_bond_type_dat_in(FILE *f,
                              const vector<AMOCH::BondType>& bond_types)
{
   unsigned i;
   int n;

   n = 1;
   for (i = 0; i < bond_types.size(); i++) {
      fprintf(f, "%d  %d  %d  %d\n", 
              n++, 
              bond_types[i].ptype1()+1, 
              bond_types[i].ptype2()+1, 
              bond_types[i].order());
   }
   fflush(f);
}

