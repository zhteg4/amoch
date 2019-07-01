// ============================================================================
// cbmc.h -- AMOCH::CBMC class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#ifndef AMOCH_CBMC_H
#define AMOCH_CBMC_H

#include "builder.h"

namespace AMOCH {

//
// CBMC - Configurational Bias Monte Carlo chain builder
// 
// See Sadanobu and Goddard, J. Chem. Phys 106, 6722-6729, 1997
//
class CBMC : public Builder {

   private:

   int _nbuild;         // number of monomers to add in build()
   int _nconfigs;       // number of MC configurations to test
   int _ndelta_steps;   // number of additional steps after configuration
   int _bonded_cutoff;  // nonbonded interactions require ? bonds between atoms
   int _nchains;        // number of Chains 
   int _bb_bondtype;    // type index of bond between monomers
   bool _bb_only;       // build only backbone particles, not side groups
   bool _rotate;        // rotate torsions to configure
   bool _only_between;  // rotate only torsions between monomers to configure
   bool _log_monomers;  // log monomers as added 
   bool _log_energies;  // log energies while configuring torsions
   bool _log_status;    // log status bar while building
   Real _inv_kT;        // 1/(kB*T)
   Real _delta_degs;    // size of additional steps (degrees)
   Real _bb_bondlen;    // length of backbone bond between monomers (Angstroms)
   std::vector<Real> _torsion_angle;     // allowed torsion values
   std::vector<Real> _torsion_energy;    // energies for allowed torsions
   std::vector<Real> _torsion_prob;      // probabilities for allowed torsions
   std::vector<int> _pattern_index;         // ChainType monomer index per Chain
   std::vector<Triple> _chain_tail_anchor;  // anchor position for each Chain
   std::vector<Triple> _chain_tail_pos;     // tail position for each Chain
   std::vector<Triple> _chain_tail_vec;     // tail vector for each Chain

   // Prepare Chain data members
   void _resize_chain_data();

   // Return the sum of nonbonded pair interactions of particles in m with
   // surrounding particles in system and exclusions
   Real _nonbond_interactions(Monomer& m, 
                              const std::vector<int>& index_map,
                              const Chain& c,
                              int ic,
                              int itail,
                              const System& sys, 
                              const std::vector<const Exclusion *>& exclusions,
                              const Energy& energy) const;

   // Return an energy for the indicated torsion angle from linear 
   // interpolation between _torsion_energy[] values
   Real _interpolate_energy(Real angle) const;

   // Set torsion angles on the backbone of the Monomer m to minimize 
   // interactions with surrounding particles
   void _rotate_torsions(Monomer& m, 
                         const std::vector<int> index_map, 
                         Chain& c, 
                         int ic, 
                         int itail, 
                         bool first, 
                         System& sys, 
                         const std::vector<const Exclusion *>& exclusions,
                         const Energy& energy,
                         RNG& rng) const;

   // Add a new bond; add to BOTH particles
   void _add_bond(System& sys,
                  Chain& ch,
                  Ivec& bond) const;

   public:

   CBMC();

   // Synthesized copy constructor, assignment, destructor

   // Setup Monomer m with parameters in pmon; return true on success or 
   // false on error after logging a message
   bool setup_monomer(Monomer& m,
                      const Param& pmon);
                              
   // Update types of particles and bonds, if any changes when forming bonds
   // between monomers m1 and m2
   void update_types(const Monomer& m1,
                     const Monomer& m2,
                     std::vector<PartType>& part_types,
                     std::vector<BondType>& bond_types);

   // Set properties from input; return true on success or false on error
   // after logging a message
   bool setup(const Param& pbuild,
              const std::vector<PartType>& part_types,
              std::vector<BondType>& bond_types);

   // Return true if chains will be built with only backbone particles, else
   // return false
   bool backbone_only() const {return _bb_only;}

   // Build chains 
   void build(System& sys,
              std::vector<Chain>& chains, 
              std::vector<Monomer>& monomers,
              std::vector<int>& monomer_selection,
              const std::vector<Chaintype>& chain_types, 
              const std::vector<const Exclusion *>& exclusions,
              const Energy& energy,
              RNG& rng);

   // Write internal state to the output FILE f
   void write(FILE *f) const;

   // Read internal state from the input FILE f
   void read(FILE *f);

};  // CBMC

}   // AMOCH

#endif  // AMOCH_CBMC_H
