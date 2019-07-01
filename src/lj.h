// ============================================================================
// lj.h -- AMOCH::LJ class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_LJ_H
#define AMOCH_LJ_H

#include "energy.h"

namespace AMOCH {

//
// Lennard-Jones interaction force field
//
class LJ : public Energy {

   private:

   std::vector<std::vector<std::vector<Real> > > _params;

   public:

   LJ(int ntypes,
      Real rcut);

   // Synthesized copy constructor, assignment, destructor

   // Set properties from input; return true on success or false on error
   bool setup(const Param& pen,
              const std::vector<PartType>& part_types);

   // Scale interaction parameters
   void scale(Real f = 1.0);

   // Calculate the pair energy between two particles of types t1, t2, 
   // separated by square distance rsq
   Real pair_energy(int t1, 
                    int t2,
                    Real rsq) const;

   // Calculate the bond energy between two particles of types t1, t2, 
   // separated by square distance rsq
   Real bond_energy(int t1, 
                    int t2,
                    Real rsq) const;

   // Calculate the hard energy between two particles of types t1, t2, 
   // separated by square distance rsq
   Real hard_energy(int t1, 
                    int t2,
                    Real rsq) const;

};  // LJ

}   // AMOCH
      
#endif  // AMOCH_LJ_H
