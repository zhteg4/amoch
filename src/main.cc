// ============================================================================
// main.cc -- AMOCH top level
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <cerrno>
#include "version.h"
#include "os.h"
#include "logfile.h"
#include "timer.h"
#include "cbmc.h"
#include "lj.h"
#include "fene.h"
#include "restart.h"
#include "output.h"
#include "slab.h"

using std::vector;
using std::string;
using AMOCH::Real;
using AMOCH::Param;
using AMOCH::Monomer;
using AMOCH::Chaintype;
using AMOCH::Chain;
using AMOCH::Triple;
using AMOCH::join_path;
using AMOCH::log_message;
using AMOCH::write_pdb;
using AMOCH::write_xyz;
using AMOCH::System;
using AMOCH::PartType;
using AMOCH::BondType;
using AMOCH::Exclusion;

// Defined in parse.cc (parse.y)
extern bool 
parse_input(Param *proot,
            const char *path);

static AMOCH::Builder *builder = NULL;
static AMOCH::Energy *energy = NULL;
static vector<const Exclusion *> exclusions(0, NULL);

// ============================================================================
// Write version and build info to output FILE f
// ============================================================================
static void
show_version(FILE *f)
{
   fprintf(f, "%s version %s\n", PROG_NAME, PROG_VERSION);
   fprintf(f, "Compiled with \"%s\"\n", PROG_CXX_CMD);
   fprintf(f, "  Linked with \"%s\"\n", PROG_LD_CMD);
   fprintf(f, "Real type: ");
#ifdef REAL_DOUBLE
   fprintf(f, "double");
#else
   fprintf(f, "float");
#endif
   fprintf(f, "\nSFMT exponent: %d\n", SFMT_MEXP);
}

// ============================================================================
// Return a pointer to a newly opened FILE, or NULL on error after logging a
// message
// ============================================================================
static FILE *
open_file(const char *path,
          const char *mode)
{
   FILE *f = fopen(path, mode);

   if (!f) 
      log_message("\nUnable to open %s: %s\n", path, strerror(errno));
   return f;
}
       
// ============================================================================
// Write files to serve as inputs to data4Lammps
// ============================================================================
static void
write_data4Lammps(const string& dir,
                  const System& sys,
                  const vector<Chain>& chains,
                  const vector<PartType>& part_types,
                  const vector<BondType>& bond_types)
{
   FILE *f;
   string path = join_path(dir, ATOMS_DATA_FILE);

   if ((f = open_file(path.c_str(), "w"))) {
      AMOCH::write_atoms_dat(f, sys, chains);
      fclose(f);
      log_message("Wrote data4Lammps input file %s\n", path.c_str());
   }
   path = join_path(dir, ATOM_TYPE_FILE);
   if ((f = open_file(path.c_str(), "w"))) {
      AMOCH::write_atom_type_dat_in(f, part_types);
      fclose(f);
      log_message("Wrote data4Lammps input file %s\n", path.c_str());
   }
   path = join_path(dir, BONDS_DATA_FILE);
   if ((f = open_file(path.c_str(), "w"))) {
      AMOCH::write_bonds_dat(f, sys, chains);
      fclose(f);
      log_message("Wrote data4Lammps input file %s\n", path.c_str());
   }
   path = join_path(dir, BOND_TYPE_FILE);
   if ((f = open_file(path.c_str(), "w"))) {
      AMOCH::write_bond_type_dat_in(f, bond_types);
      fclose(f);
      log_message("Wrote data4Lammps input file %s\n", path.c_str());
   }
}

// ============================================================================
// Prepare to exit
// ============================================================================
static void
cleanup()
{
   if ((builder))
      delete builder;
   if ((energy))
      delete energy;
   for (unsigned i = 0; i < exclusions.size(); i++)
      delete exclusions[i];
   AMOCH::close_logfile();
}

// ============================================================================
// Return EXIT_SUCCESS or EXIT_FAILURE
// ============================================================================
int
main(int argc,
     char **argv)
{
   // -------------------------------------------------
   // Initialize
   // -------------------------------------------------
   AMOCH::Timer timer;
   FILE *f;
   int i, j, k, N;

   timer.start();
   atexit(cleanup);

   // -------------------------------------------------
   // Command line
   // -------------------------------------------------
   const char usage[] = 
      "Usage: %s [flags] infile\n"
      "Flags:\n"
      "   -s   Setup only\n"
      "   -v   Show version and build info, then quit\n";
#define SHOW_USAGE  fprintf(stderr, usage, PROG_NAME)
   char *infile = NULL;
   char *p = NULL;
   bool setup_only = false;

   while (--argc > 0) {
      p = *++argv;
      if ('-' == *p) {
         ++p;
         if (*(p+1)) {
            SHOW_USAGE;
            fprintf(stderr, "\nUnknown flag: \"%s\"\n", p);
            exit(EXIT_FAILURE);
         }
         switch (*p) {
            case 's':
               setup_only = true;
               break;
            case 'v':
               show_version(stdout);
               return EXIT_SUCCESS;
               break;
            default:
               SHOW_USAGE;
               fprintf(stderr, "\nUnknown flag: '%c'\n", *p);
               exit(EXIT_FAILURE);
               break;
         }
      }  // end flag
      else {
         if ((infile)) {
            SHOW_USAGE;
            fprintf(stderr, "\nMultiple input files: %s, %s\n", infile, p);
            exit(EXIT_FAILURE);
         }
         infile = p;
      }
   }
   if (!infile) {
      SHOW_USAGE;
      fprintf(stderr, "\nMissing input file\n");
      exit(EXIT_FAILURE);
   }

   // -------------------------------------------------
   // Parse input
   // -------------------------------------------------
   Param proot("root");
   const Param *pin, *pin2;
   string includedir("");
   string path;

   if (!parse_input(&proot, infile))  // message written to stderr
      exit(EXIT_FAILURE);

   if ((pin = proot.find("includedir")))
      includedir = pin->value();

   pin = proot.find("include");
   while ((pin)) {
      path = pin->value();
      if (includedir.size() > 0) 
         path = join_path(includedir, pin->value());
      if (!parse_input(&proot, path.c_str()))
         exit(EXIT_FAILURE);
      pin = pin->find_next("include");
   }

   // -------------------------------------------------
   // Logging -- use log_message() after this
   // -------------------------------------------------
   AMOCH::set_logfile(stdout);

   if ((pin = proot.find("write_log")))
      AMOCH::open_logfile(pin->value().c_str());
   log_message("\n\n");
   show_version(AMOCH::get_logfile());
   log_message("\n\n");

   // -------------------------------------------------
   // Restart -- initial
   // -------------------------------------------------
   string restartdir("");

   if ((pin = proot.find("restart")))
      restartdir = pin->value();

   // -------------------------------------------------
   // RNG 
   // -------------------------------------------------
   uint32_t rng_seed = time(NULL);

   if ((pin = proot.find("rng_seed")))
      rng_seed = (uint32_t)atoi(pin->value().c_str());
   else {
      if (restartdir.size() > 0) {
         log_message("\nWARNING: restarting simulation without RNG seed!\n");
         exit(EXIT_FAILURE);
      }
   }
   log_message("\nRNG seed: %u\n\n", rng_seed);

   AMOCH::RNG rng(rng_seed);

   // -------------------------------------------------
   // Builder -- initial
   // -------------------------------------------------
   const Param *pbuild = proot.find("builder");

   if ((pbuild)) {
      if (pbuild->value() == "cbmc")
         builder = new AMOCH::CBMC;
      // XXX more Builders here
      else {
         log_message("\nUnknown builder: %s\n", pbuild->value().c_str());
         exit(EXIT_FAILURE);
      }
   }
   else {
      log_message("\nNo builder defined\n");
      exit(EXIT_FAILURE);
   }

   // -------------------------------------------------
   // Monomers
   // -------------------------------------------------
   vector<Monomer> monomers(0);
   vector<PartType> part_types(0);
   vector<BondType> bond_types(0);
   const Param *pmon = NULL;
   int max_bonds_per_particle = 0;
   int mmb;

   OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obInfo);
   pmon = proot.find("monomer");
   while ((pmon)) {
      monomers.push_back(Monomer(pmon->value()));
      Monomer& m = monomers.back();

      if ((pin = pmon->find("file"))) {
         path = pin->value();
         if (!AMOCH::check_file(path) && (includedir.size() > 0)) 
            path = join_path(includedir, pin->value());
         if (!m.read(path.c_str()))  // message logged
            exit(EXIT_FAILURE);
      }
      else {
         log_message("\nNo file specified for monomer %s\n", m.name().c_str());
         exit(EXIT_FAILURE);
      }

      if (!builder->setup_monomer(m, *pmon))  // message logged
         exit(EXIT_FAILURE);

      if (!m.setup(*pmon, part_types))  // message logged
         exit(EXIT_FAILURE);
      m.find_bonds(bond_types);  // within the molecule m
      m.adjust_internal_coords(part_types);

      mmb = m.max_bonds_per_atom();
      if (mmb > max_bonds_per_particle)
         max_bonds_per_particle = mmb;

      pmon = pmon->find_next("monomer");
      log_message("\n");
   }
   if (0 == monomers.size()) {
      log_message("\nNo monomer specified\n");
      exit(EXIT_FAILURE);
   }
   log_message("After reading all monomers, found %u particle types:\n", 
               part_types.size());
   for (i = 0; i < (int)part_types.size(); i++) {
      log_message("   %d. %s with %d bonds", i+1,part_types[i].symbol().c_str(),
                  part_types[i].nbonds());
      if (part_types[i].hbond())
         log_message(" (H bond)");
      if (part_types[i].resonant())
         log_message(" (resonant)");
      log_message("\n");
   }
   log_message("\nAfter reading all monomers, found %u bond types:\n", 
               bond_types.size());
   for (i = 0; i < (int)bond_types.size(); i++) {
      log_message("   %d. %d - %d, order %d\n", i+1, bond_types[i].ptype1()+1,
                  bond_types[i].ptype2()+1, bond_types[i].order());
   }

   // -------------------------------------------------
   // Builder
   // -------------------------------------------------
   if (!builder->setup(*pbuild, part_types, bond_types))  // message logged
      exit(EXIT_FAILURE);
   log_message("\n");

   // -------------------------------------------------
   // Chain types
   // -------------------------------------------------
   vector<Chaintype> chain_types(0);
   vector<Real> chain_type_weights(0);
   Real wt;

   pin = proot.find("chain_type");
   while ((pin)) {
      Chaintype ct(pin->value());

      if (!ct.setup(*pin, monomers))  // message logged
         exit(EXIT_FAILURE);
      chain_types.push_back(ct);

      if (ct.select_by_weight()) {
         for (i = 0; i < ct.nmonomers(); i++) {
            for (j = 0; j < ct.nmonomers(); j++) {
               builder->update_types(monomers[ct.monomer_index(i)], 
                                     monomers[ct.monomer_index(j)],
                                     part_types, bond_types);
            }
         }
      }
      else {
         for (i = 0; i < ct.nmonomers(); i++) {
            for (j = i; j < ct.nmonomers(); j++) {
               builder->update_types(monomers[ct.monomer_index(i)], 
                                     monomers[ct.monomer_index(j)],
                                     part_types, bond_types);
            }
         }
      }

      if ((pin2 = pin->find("weight"))) {
         wt = atof(pin2->value().c_str());
         chain_type_weights.push_back(wt);
         log_message("   chain type weight: %g\n", wt);
      }
      pin = pin->find_next("chain_type");
      log_message("\n");
   }
   N = chain_types.size();
   if (0 == N) {
      log_message("\nNo chain types specified\n");
      exit(EXIT_FAILURE);
   }
   if (0 == chain_type_weights.size()) {
      log_message("No chain_type weights specified; ");
      log_message("all chain_types have equal weight\n");
      wt = 1.0/(Real)N;
      for (i = 0; i < N; i++)
         chain_type_weights.push_back(wt);
   }
   wt = 0.0;
   for (i = 0; i < N; i++)
      wt += chain_type_weights[i];
   if (fabs(1.0-wt) > REAL_SMALL) {
      log_message("\nSum of chain_type weights != 1\n");
      exit(EXIT_FAILURE);
   }

   // -------------------------------------------------
   // Chains
   // -------------------------------------------------
   vector<Chain> chains(0);
   int type, Nchains, Nmonomers;
   Real stddev = 0.0;
   bool have_lengths = false;

   pin = proot.find("chains");
   if (!pin) {
      log_message("\nChain length distribution not specified\n");
      exit(EXIT_FAILURE);
   }
   Nchains = atoi(pin->value().c_str());
   log_message("Build %d chains:\n", Nchains);
   chains.reserve(Nchains);

   pin2 = pin->find("monomers");
   if (!pin2) {
      log_message("\nMonomers per chain (polymerization) not specified\n");
      exit(EXIT_FAILURE);
   }
   Nmonomers = atoi(pin2->value().c_str());
   log_message("   Mean polymerization: %d monomers per chain\n", Nmonomers);
   if ((pin2 = pin->find("stddev"))) {
      stddev = atof(pin2->value().c_str());
      log_message("   Standard deviation in chain length: %g\n", stddev);
   }
   else 
      log_message("   All chains have the same length\n");
   rng.set_normal((Real)Nmonomers, stddev);

   if ((pin2 = pin->find("lengths"))) {
      log_message("   Chain lengths specified\n");
      have_lengths = true;
      N = pin2->nvalues();
      if (N != Nchains) {
         log_message("\nChain lengths specified for %d chains; ", N);
         log_message("expecting %d chains\n", Nchains);
         exit(EXIT_FAILURE);
      }
   }

   int chain_length_sum = 0;
   int total_particles = 0;
   int nm, np;
   vector<int> chain_type_selection(chain_types.size(), 0);
   Real chmass;
   Real total_mass = 0.0;

   for (i = 0; i < Nchains; i++) {
      type = rng.select_by_weight(chain_type_weights);
      chain_type_selection[type]++;
      if (have_lengths) 
         Nmonomers = atoi(pin2->value(i).c_str());
      else
         Nmonomers = rng.yield_normal();
      chain_length_sum += Nmonomers;

      chmass = 0.0;
      np = 0;
      N = chain_types[type].nmonomers();
      nm = (int)ceil((Real)Nmonomers/(Real)N);
      for (j = 0; j < N; j++) {
         if (chain_types[type].select_by_weight())
            nm = (int)ceil(chain_types[type].monomer_weight(j) * Nmonomers);
         // nm is an estimate of how many times monomer j of ChainType type 
         // will appear in Chain i (out of Nmonomers monomers)

         if (builder->backbone_only()) {
            chmass += nm*monomers[chain_types[type].monomer_index(j)].bbmass();
            np += nm * monomers[chain_types[type].monomer_index(j)].nbackbone();
         }
         else {
            chmass += nm * monomers[chain_types[type].monomer_index(j)].mass();
            np += nm * monomers[chain_types[type].monomer_index(j)].natoms();
         }
      }
      total_particles += np;
      total_mass += chmass;
      chains.push_back(Chain(type, Nmonomers, np));
   }
   log_message("   Actual mean polymerization: %g monomers per chain\n", 
               (Real)chain_length_sum / (Real)Nchains);
   log_message("   Chain type distribution:\n");
   for (i = 0; i < (int)chain_types.size(); i++) {
      log_message("      %s: %g\n", chain_types[i].name().c_str(), 
                  (Real)chain_type_selection[i] / (Real)Nchains);
   }
   log_message("\n");
   if (0 == total_particles) {
      log_message("\nNo particles in simulation\n");
      exit(EXIT_FAILURE);
   }

   // -------------------------------------------------
   // Energy expression
   // -------------------------------------------------
   Real rcut = 0.0;

   pin = proot.find("energy_cutoff");
   if (!pin) {
      log_message("\nEnergy cutoff not specified\n");
      exit(EXIT_FAILURE);
   }
   rcut = atof(pin->value().c_str());

   pin = proot.find("energy");
   if (!pin) {
      log_message("\nEnergy expression not specified\n");
      exit(EXIT_FAILURE);
   }
   if (pin->value() == "LJ")
      energy = new AMOCH::LJ(part_types.size(), rcut);
   else if (pin->value() == "FENE")
      energy = new AMOCH::FENE(part_types.size(), rcut);
   // XXX more Energy types here
   else {
      log_message("\nUnknown energy expression: %s\n", pin->value().c_str());
      exit(EXIT_FAILURE);
   }
   log_message("%s energy expression:\n", pin->value().c_str());
   log_message("   pair interaction cutoff: %g Angstroms\n", rcut);
   if (!energy->setup(*pin, part_types))  // message logged
      exit(EXIT_FAILURE);

   Real escale = 1.0;

   if ((pin = proot.find("energy_scale"))) {
      escale = atof(pin->value().c_str());
      log_message("   Scale interaction parameters by %g\n", escale);
   }
   energy->scale(escale);

   // -------------------------------------------------
   // System
   // -------------------------------------------------
   Triple min, max;
   const Real density_scale = 1.6605;
   const Param *pexc = NULL;
   Real density, volume;
   Exclusion *ex = NULL;

   log_message("\n%d particles in simulation\n", total_particles);
   log_message("Simulation volume:\n");
   min.zero();
   if ((pin = proot.find("density"))) {
      density = atof(pin->value().c_str());  // g / cm^3
      volume = (total_mass / density) * density_scale;  // A^3
      max.x = cbrt(volume);
      max.y = max.x;
      max.z = max.x;
      log_message("   Requested density: %g g/cm^3\n", density);
      log_message("   Edge length: %g Angstroms\n", max.x);
   }
   else if ((pin = proot.find("system"))) {
      if (3 != pin->nvalues()) {
         log_message("\nExpected [x, y, z] system size\n");
         exit(EXIT_FAILURE);
      }
      max.x = atof(pin->value(0).c_str());
      max.y = atof(pin->value(1).c_str());
      max.z = atof(pin->value(2).c_str());
      volume = max.x * max.y * max.z;  // A^3
      log_message("   Dimensions %g x %g x %g Angstroms\n", max.x, max.y,max.z);

      pexc = proot.find("exclude");
      while ((pexc)) {
         ex = NULL;
         if (pexc->value() == "slab")
            ex = new AMOCH::Slab;
         // XXX more excluded volumes here
         else {
            log_message("\nIgnoring unknown excluded volume: %s\n", 
                        pexc->value().c_str());
         }
         if ((ex)) {
            if (!ex->setup(*pexc))  // message logged
               exit(EXIT_FAILURE);
            volume -= ex->volume();
            exclusions.push_back(ex);
         }
         pexc = pexc->find_next("exclude");
      }
   }
   else {
      log_message("\nNeither density nor system dimensions specified\n");
      exit(EXIT_FAILURE);
   }
   density = (total_mass / volume) * density_scale;  // g / cm^3
   log_message("   Estimated density: %g g/cm^3\n", density);

   const Real cell_size = sqrt(energy->pair_range_sq()) + RANGE_EPS;
   System sys(min, max, cell_size, total_particles, 
                     max_bonds_per_particle);

   log_message("   %d cells\n", sys.ncells());

   // -------------------------------------------------
   // Restart
   // -------------------------------------------------
   const int pathlen = 128;
   char cpath[pathlen];

   if (restartdir.size() > 0) {
      bool ret;

      log_message("\nRestart simulation from %s:\n", restartdir.c_str());

      // RNG
      path = join_path(restartdir, RNG_DATA_FILE);
      if (!(f = open_file(path.c_str(), "r")))  // message logged
         exit(EXIT_FAILURE);
      rng.read(f);
      fclose(f);
      log_message("   RNG state\n");

      // Chains
      for (i = 0; i < Nchains; i++) {
         snprintf(cpath, pathlen, CHAIN_DATA_FILE, i);
         path = join_path(restartdir, cpath);
         if (!(f = open_file(path.c_str(), "r")))  // message logged
            exit(EXIT_FAILURE);
         chains[i].read(f);
         fclose(f);
      }
      log_message("   Chain states\n");

      // LAMMPS dump - particles
      path = join_path(restartdir, LAMMPS_DUMP_FILE);
      if (!(f = open_file(path.c_str(), "r")))  // message logged
         exit(EXIT_FAILURE);
      ret = AMOCH::read_lammps_dump(f, sys, chains);
      fclose(f);
      if (!ret)  // message logged
         exit(EXIT_FAILURE);
      log_message("   LAMMPS particle dump\n");

      // Bonds
      path = join_path(restartdir, BONDS_DATA_FILE);
      if (!(f = open_file(path.c_str(), "r")))  // message logged
         exit(EXIT_FAILURE);
      ret = AMOCH::read_bonds(f, sys, chains);
      fclose(f);
      if (!ret)  // message logged
         exit(EXIT_FAILURE);
      log_message("   Bonds\n");

      // Builder
      path = join_path(restartdir, BUILDER_DATA_FILE);
      if (!(f = open_file(path.c_str(), "r")))  // message logged
         exit(EXIT_FAILURE);
      builder->read(f);
      fclose(f);
      log_message("   Builder\n\n");
   }  // end restart

   // -------------------------------------------------
   // End setup
   // -------------------------------------------------
   Real setup_time = timer.elapsed();

   log_message("\n\n");
   log_message("Setup complete in %.6f seconds\n", setup_time);
   log_message("\n\n");
   if (setup_only)
      return EXIT_SUCCESS;

   // -------------------------------------------------
   // Build chains
   // -------------------------------------------------
   vector<int> monomer_selection(monomers.size(), 0);

   builder->build(sys, chains, monomers, monomer_selection, chain_types, 
                  exclusions, *energy, rng);

   Real build_time = timer.elapsed() - setup_time;

   log_message("\n\n");
   log_message("Build complete in %.6f seconds\n", build_time);
   log_message("\n\n");

   // -------------------------------------------------
   // Output -- monomer selections
   // -------------------------------------------------
   int n = 0;
   for (i = 0; i < (int)monomer_selection.size(); i++)
      n += monomer_selection[i];
   log_message("Monomer selections:\n");
   for (i = 0; i < (int)monomer_selection.size(); i++) {
      log_message("   %s: %d (%.3f)\n", monomers[i].name().c_str(), 
                  monomer_selection[i], (Real)monomer_selection[i]/(Real)n);
   }
   log_message("\n");

   // -------------------------------------------------
   // Output -- molecular structures
   // -------------------------------------------------
   const string defcomm("Chains built using "PROG_NAME" version "PROG_VERSION);
   string comment;
   bool connect;

   if ((pin = proot.find("write_wrapped_pdb"))) {
      path = pin->value();
      if ((f = open_file(path.c_str(), "w"))) {
         comment = ((pin2 = pin->find("comment"))) ? pin2->value() : defcomm;
         connect = true;
         if ((pin2 = pin->find("exclude"))) {
            if (pin2->value() == "connect")
               connect = false;
            else  {
               log_message("\nIgnoring unknown \"exclude\" value: %s\n", 
                           pin2->value().c_str());
            }
         }
         write_pdb(f, comment.c_str(), sys, chains, part_types, true, connect);
         fclose(f);
         log_message("Wrote wrapped PDB structure file %s\n", path.c_str());
      }
   }
   if ((pin = proot.find("write_unwrapped_pdb"))) {
      path = pin->value();
      if ((f = open_file(path.c_str(), "w"))) {
         comment = ((pin2 = pin->find("comment"))) ? pin2->value() : defcomm;
         connect = true;
         if ((pin2 = pin->find("exclude"))) {
            if (pin2->value() == "connect")
               connect = false;
            else  {
               log_message("\nIgnoring unknown \"exclude\" value: %s\n", 
                           pin2->value().c_str());
            }
         }
         write_pdb(f, comment.c_str(), sys, chains, part_types, false, connect);
         fclose(f);
         log_message("Wrote unwrapped PDB structure file %s\n", path.c_str());
      }
   }

   if ((pin = proot.find("write_wrapped_xyz"))) {
      path = pin->value();
      if ((f = open_file(path.c_str(), "w"))) {
         comment = ((pin2 = pin->find("comment"))) ? pin2->value() : defcomm;
         write_xyz(f, comment.c_str(), sys, chains, part_types, true);
         fclose(f);
         log_message("Wrote wrapped XYZ structure file %s\n", path.c_str());
      }
   }
   if ((pin = proot.find("write_unwrapped_xyz"))) {
      path = pin->value();
      if ((f = open_file(path.c_str(), "w"))) {
         comment = ((pin2 = pin->find("comment"))) ? pin2->value() : defcomm;
         write_xyz(f, comment.c_str(), sys, chains, part_types, false);
         fclose(f);
         log_message("Wrote unwrapped XYZ structure file %s\n", path.c_str());
      }
   }

   // -------------------------------------------------
   // Output -- chain statistics
   // -------------------------------------------------
   const Param *plen = proot.find("write_chain_length");
   const Param *phist = proot.find("write_chain_length_histo");
   const Param *pmsid = proot.find("write_msid");
   const Param *pmsid_atom = proot.find("write_msid_atom");
   int cell, slot;
   int max_monomers = 0;
   int max_particles = 0;
   Real max_length = 0.0;
   Real bin_size = 5.0;  // Angstroms
   Real len;
   Triple tail_pos;

   if ((plen) || (phist) || (pmsid) || (pmsid_atom)) {
      log_message("\nGenerating chain statistics\n\n");
      for (i = 0; i < Nchains; i++) {
         for (j = 0; j < chains[i].nparticles(); j++) {
            cell = chains[i].cell(j);
            slot = chains[i].slot(j);
            chains[i].record_particle(j, sys.cells[cell].tpos[slot]);
         }
         j = chains[i].tail_index;
         cell = chains[i].cell(j);
         slot = chains[i].slot(j);
         tail_pos = sys.cells[cell].tpos[slot];
         chains[i].finalize(tail_pos);
         if (chains[i].max_monomers() > max_monomers)
            max_monomers = chains[i].max_monomers();
         if (chains[i].nparticles() > max_particles)
            max_particles = chains[i].nparticles();
      }
   }

   vector<Real> length(max_monomers, 0.0);
   vector<int> count(max_monomers, 0);

   // Chain length stats
   if ((plen)) {
      for (i = 0; i < Nchains; i++) {
         for (j = 0; j < chains[i].max_monomers(); j++) {
            len = chains[i].length(j);
            length[j] += len;
            count[j]++;
            if (len > max_length)
               max_length = len;
         }
      }
      path = plen->value();
      if ((f = open_file(path.c_str(), "w"))) {
         for (j = 0; j < max_monomers; j++)
            fprintf(f, "%5d%8.2f\n", j, length[j]/(Real)count[j]);
         fclose(f);
         log_message("Wrote chain length file %s\n", path.c_str());
      }
   }

   // Chain length histogram
   if ((phist)) {
      if ((pin2 = phist->find("bin_size")))
         bin_size = atof(pin2->value().c_str());

      int nbins = (int)ceil(max_length/bin_size);
      vector<int> lenhisto(nbins, 0);

      for (i = 0; i < Nchains; i++) {
         len = chains[i].length(chains[i].max_monomers() - 1);
         lenhisto[(int)(len/bin_size)]++;
      }
      path = phist->value();
      if ((f = open_file(path.c_str(), "w"))) {
         for (i = 0; i < nbins; i++)
            fprintf(f, "%g %d\n", (i+1)*bin_size, lenhisto[i]);
         fclose(f);
         log_message("Wrote chain length histogram file %s ", path.c_str());
         log_message("with bin size %g Ang\n", bin_size);
      }
   }

   // Mean square internal displacements between monomers
   if ((pmsid)) {
      for (i = 0; i < max_monomers; i++) {
         length[i] = 0.0;
         count[i] = 0;
      }
      for (i = 0; i < Nchains; i++) {
         for (j = 0; j < chains[i].max_monomers()-1; j++) {
            for (k = j+1; k < chains[i].max_monomers(); k++) {
               length[k-j] += (chains[i].com(k) - chains[i].com(j)).sqlength();
               count[k-j]++;
            }
         }
      }
      path = pmsid->value();
      if ((f = open_file(path.c_str(), "w"))) {
         Real g;

         for (i = 1; i < max_monomers; i++) {
            g = length[i]/(Real)count[i];
            fprintf(f, "%5d%8.2f%8.2f\n", i, g, g/(Real)i);
         }
         fclose(f);
         log_message("Wrote monomer msid file %s\n", path.c_str());
      }
   }

   // Mean square internal displacements between atoms
   if ((pmsid_atom)) {
      vector<Real> alength(max_particles, 0.0);
      vector<int> acount(max_particles, 0);
      Triple pj, pk;

      for (i = 0; i < Nchains; i++) {
         for (j = 0; j < chains[i].nparticles()-1; j++) {
            cell = chains[i].cell(j);
            slot = chains[i].slot(j);
            pj = sys.cells[cell].tpos[slot];
            for (k = j+1; k < chains[i].nparticles(); k++) {
               cell = chains[i].cell(k);
               slot = chains[i].slot(k);
               pk = sys.cells[cell].tpos[slot];
               alength[k-j] += (pk - pj).sqlength();
               acount[k-j]++;
            }
         }
      }
      path = pmsid_atom->value();
      if ((f = open_file(path.c_str(), "w"))) {
         Real g;

         for (i = 1; i < max_particles; i++) {
            g = alength[i]/(Real)acount[i];
            fprintf(f, "%5d%8.2f%8.2f\n", i, g, g/(Real)i);
         }
         fclose(f);
         log_message("Wrote atomic msid file %s\n", path.c_str());
      }
   }

   // Torsion histogram
   if ((pin = proot.find("write_torsion_histo"))) {
      path = pin->value();
      if ((f = open_file(path.c_str(), "w"))) {
         vector<int> torsion_count(360, 0);

         for (i = 0; i < Nchains; i++) {
            for (j = 0; j < 360; j++)
               torsion_count[j] += chains[i].torsion_count(j);
         }
         for (j = 0; j < 360; j++) {
            fprintf(f, "%d  %d\n", j, torsion_count[j]);
         }
         fclose(f);
         log_message("Wrote torsion histogram file %s\n", path.c_str());
      }
   }

   // -------------------------------------------------
   // Output -- data4Lammps
   // -------------------------------------------------
   bool write = true;

   if ((pin = proot.find("write_data4Lammps"))) {
      string dir = pin->value();

      if (!AMOCH::check_dir(dir)) {
         if (!AMOCH::create_dir(dir))
            write = false;
      }
      if (write)
         write_data4Lammps(dir, sys, chains, part_types, bond_types);
   }

   // -------------------------------------------------
   // Output - restart
   // -------------------------------------------------
   if ((pin = proot.find("write_restart"))) {
      restartdir = pin->value();
      log_message("\nWrite restart files to %s\n\n", restartdir.c_str());
      write = true;

      if (!AMOCH::check_dir(restartdir)) {
         if (!AMOCH::create_dir(restartdir))
            write = false;
      }
      if (write) {
         write_data4Lammps(restartdir, sys, chains, part_types, bond_types);

         for (i = 0; i < Nchains; i++) {
            snprintf(cpath, pathlen, CHAIN_DATA_FILE, i);
            path = join_path(restartdir, cpath);
            if ((f = open_file(path.c_str(), "w"))) {
               chains[i].write(f);
               fclose(f);
            }
         }
         log_message("Wrote chain state files to %s\n", restartdir.c_str());

         path = join_path(restartdir, RNG_DATA_FILE);
         if ((f = open_file(path.c_str(), "w"))) {
            rng.write(f);
            fclose(f);
         }
         log_message("Wrote RNG state to %s\n", restartdir.c_str());

         path = join_path(restartdir, BUILDER_DATA_FILE);
         if ((f = open_file(path.c_str(), "w"))) {
            builder->write(f);
            fclose(f);
         }
         log_message("Wrote builder state to %s\n", restartdir.c_str());
      }
   }

   // -------------------------------------------------
   // Output - miscellaneous
   // -------------------------------------------------
   if ((pin = proot.find("write_cell"))) {
      if ((f = open_file(pin->value().c_str(), "w"))) {
         fprintf(f, "%g  %g  %g\n", max.x, max.y, max.z);
         fclose(f);
         log_message("Wrote cell file %s\n", pin->value().c_str());
      }
   }

   Real output_time = timer.elapsed() - (setup_time + build_time);

   // -------------------------------------------------
   // Finish 
   // -------------------------------------------------
   log_message("\nTotal simulation time: %.6f seconds\n", timer.elapsed());
   log_message("   Setup time:  %.6f seconds\n", setup_time);
   log_message("   Build time:  %.6f seconds\n", build_time);
   log_message("      Cell expansions: %d\n", AMOCH::Cell::nresize);
   log_message("   Output time: %.6f seconds\n\n", output_time);
   return EXIT_SUCCESS;
}
