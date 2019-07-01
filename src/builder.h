// ============================================================================
// builder.h -- AMOCH::Builder class
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#ifndef AMOCH_BUILDER_H
#define AMOCH_BUILDER_H

#include "system.h"
#include "chaintype.h"
#include "exclusion.h"
#include "energy.h"

namespace AMOCH {

//
// Base class for specific builder classes
//
class Builder {

   public:

   Builder() {}

   // Synthesized copy constructor, assignment

   virtual ~Builder() {}

   // Setup Monomer m with parameters in pmon; return true on success or 
   // false on error after logging a message
   virtual bool setup_monomer(Monomer& m,
                              const Param& pmon) = 0;
                              
   // Update types of particles and bonds, if any changes when forming bonds
   // between monomers m1 and m2
   virtual void update_types(const Monomer& m1,
                             const Monomer& m2,
                             std::vector<PartType>& part_types,
                             std::vector<BondType>& bond_types) = 0;

   // Set properties from input; return true on success or false on error
   // after logging a message
   virtual bool setup(const Param& pbuild,
                      const std::vector<PartType>& part_types,
                      std::vector<BondType>& bond_types) = 0;

   // Return true if chains will be built with only backbone particles, else
   // return false
   virtual bool backbone_only() const = 0;

   // Build chains 
   virtual void build(System& sys,
                      std::vector<Chain>& chains, 
                      std::vector<Monomer>& monomers,
                      std::vector<int>& monomer_selection,
                      const std::vector<Chaintype>& chain_types, 
                      const std::vector<const Exclusion *>& exclusions,
                      const Energy& energy,
                      RNG& rng) = 0;

   // Write internal state to the output FILE f
   virtual void write(FILE *f) const = 0;

   // Read internal state from the input FILE f
   virtual void read(FILE *f) = 0;

};  // Builder

}   // AMOCH

#endif  // AMOCH_BUILDER_H
