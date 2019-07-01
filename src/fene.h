// ============================================================================
// fene.h -- AMOCH::FENE class
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_FENE_H
#define AMOCH_FENE_H

#include "energy.h"

namespace AMOCH {

//
// Finite Extensible Nonlinear Elastic interaction force field
//
// See J. Chem. Phys. v119, 12718 (2003)
//
class FENE : public Energy {

   private:

   std::vector<std::vector<std::vector<Real> > > _params;

   public:

   FENE(int ntypes,
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

};  // FENE

}   // AMOCH
      
#endif  // AMOCH_FENE_H
