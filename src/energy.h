// ============================================================================
// energy.h -- AMOCH::Energy class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_ENERGY_H
#define AMOCH_ENERGY_H

#include "config.h"
#include "param.h"
#include "types.h"

namespace AMOCH {

//
// Particle interaction energy expression
//
// Abstract base class for LJ, FENE
//
class Energy {

   protected:

   int _ntypes;            // number of particle species
   int _nparams_in;        // number of force field parameters per species
   Real _pair_range_sq;    // square pair cutoff  (Ang^2)
   std::vector<std::vector<Real> > _params_in;   // input parameters

   public:

   Energy(int ntypes,
          int nparams,
          Real rcut) : 
      _ntypes(ntypes),
      _nparams_in(nparams),
      _pair_range_sq(rcut*rcut),
      _params_in(ntypes, std::vector<Real>(nparams))
   {}

   // Synthesized copy constructor, assignment

   virtual ~Energy() {}

   // Set properties from input; return true on success or false on error
   virtual bool setup(const Param& pen,
                      const std::vector<PartType>& part_types) = 0;

   // Scale interaction parameters
   virtual void scale(Real f = 1.0) = 0;

   // Return the square pair interaction range
   Real pair_range_sq() const {return _pair_range_sq;}

   // Calculate the pair energy between two particles of types t1, t2, 
   // separated by square distance rsq
   virtual Real pair_energy(int t1, 
                            int t2,
                            Real rsq) const = 0;

   // Calculate the bond energy between two particles of types t1, t2, 
   // separated by square distance rsq
   virtual Real bond_energy(int t1, 
                            int t2,
                            Real rsq) const = 0;

   // Calculate the hard energy between two particles of types t1, t2, 
   // separated by square distance rsq
   virtual Real hard_energy(int t1, 
                            int t2,
                            Real rsq) const = 0;

};  // Energy

}   // AMOCH
      
#endif  // AMOCH_ENERGY_H
