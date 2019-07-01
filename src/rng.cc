// ============================================================================
// rng.cc -- AMOCH::RNG methods
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#include <cmath>
#include "rng.h"

using AMOCH::RNG;
using AMOCH::Real;

// ============================================================================
// Constructor
// ============================================================================
RNG::RNG(uint32_t seed) : 
   _mean(0.0),
   _stddev(0.0),
   _cached_val(0.0),
   _use_cached(false)
{
   sfmt_init_gen_rand(&_state, seed);
}

// ============================================================================
// Set mean and standard deviation for yield_normal() calls
// ============================================================================
void 
RNG::set_normal(Real mean,
                Real stddev)
{
   _mean = mean;
   _stddev = stddev;
   _use_cached = false;
}

// ============================================================================
// Return a pseudorandom Real from a normal distribution with mean and
// standard deviation passed to set_normal()
// ============================================================================
Real 
RNG::yield_normal()
{
   Real retval;

   if (_use_cached) {
      retval = _cached_val;
      _use_cached = false;
   }
   else {
      Real f1, f2, z;

      // Polar Box-Muller
      do {
         f1 = 2.0*yield() - 1.0;
         f2 = 2.0*yield() - 1.0;
          z = f1*f1 + f2*f2;
      } while (z >= 1.0);
      z = sqrt(-2.0*log(z)/z);
      retval = _mean + f1*z*_stddev;
      _cached_val = _mean + f2*z*_stddev;
      _use_cached = true;
   }
   return retval;
}

// ============================================================================
// Return an index in [0, nweights) chosen according to weights
// ============================================================================
int 
RNG::select_by_weight(const std::vector<Real>& weights)
{
   const Real target = yield();
   Real p = 0.0;
   unsigned i;

   for (i = 0; i < weights.size()-1; i++) {
      if ((target >= p) && (target < p+weights[i]))
         break;
      p += weights[i];
   }
   return i;
}

// ============================================================================
// Write state to the output FILE f
// ============================================================================
void
RNG::write(FILE *f)
{
   fwrite((void *)&_state.state[0], sizeof(w128_t), SFMT_N, f);
   fwrite((void *)&_state.idx, sizeof(int), 1, f);
   fwrite((void *)&_mean, sizeof(Real), 1, f);
   fwrite((void *)&_stddev, sizeof(Real), 1, f);
   fwrite((void *)&_cached_val, sizeof(Real), 1, f);
   fwrite((void *)&_use_cached, sizeof(int), 1, f);
}

// ============================================================================
// Read state from the input FILE f
// ============================================================================
void
RNG::read(FILE *f)
{
   fread((void *)&_state.state[0], sizeof(w128_t), SFMT_N, f);
   fread((void *)&_state.idx, sizeof(int), 1, f);
   fread((void *)&_mean, sizeof(Real), 1, f);
   fread((void *)&_stddev, sizeof(Real), 1, f);
   fread((void *)&_cached_val, sizeof(Real), 1, f);
   fread((void *)&_use_cached, sizeof(int), 1, f);
}

