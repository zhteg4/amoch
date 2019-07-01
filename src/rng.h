// ============================================================================
// rng.h -- AMOCH::RNG class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#ifndef AMOCH_RNG_H  
#define AMOCH_RNG_H

#include <vector>
#include <cstdio>
#include <stdint.h>
#undef HAVE_ALTIVEC
#undef HAVE_SSE2
#include "SFMT/SFMT.h"

namespace AMOCH {

//
// Pseudorandom number generator (wrapper around SIMD fast Mersenne Twister)
//
class RNG {

   private:

   sfmt_t _state;
   Real _mean;
   Real _stddev;
   Real _cached_val;
   bool _use_cached;

   public:

   explicit RNG(uint32_t seed = 0);

   // Synthesized copy constructor, assignment, destructor

   // Return a pseudorandom Real uniformly distributed on [0,1).
   Real yield() {return sfmt_genrand_real2(&_state);}
   
   // Set mean and standard deviation for yield_normal() calls
   void set_normal(Real mean,
                   Real stddev);

   // Return a pseudorandom Real from a normal distribution with mean and
   // standard deviation passed to set_normal()
   Real yield_normal();

   // Return an index in [0, nweights) chosen according to weights
   int select_by_weight(const std::vector<Real>& weights);

   // Write state to the output FILE f
   void write(FILE *f);

   // Read state from the input FILE f
   void read(FILE *f);

};  // RNG

}   // AMOCH

#endif  // AMOCH_RNG_H
