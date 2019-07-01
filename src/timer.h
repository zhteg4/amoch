// ============================================================================
// timer.h -- AMOCH::Timer class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_TIMER_H
#define AMOCH_TIMER_H

#include <sys/time.h>
#include "config.h"

namespace AMOCH {

//
// Measure time to microseconds
//
class Timer {

   private:

   Real _t0;
   struct timeval _tv;

   // Return the current time to microseconds
   Real _now()
   {
      gettimeofday(&_tv, NULL);
      return (Real)_tv.tv_sec + 1e-6*(Real)_tv.tv_usec;
   }

   public:

   Timer() : _t0(0.0) {}

   // Synthesized copy constructor, assignment, destructor

   // Mark start time
   void start() {_t0 = _now();}

   // Return the elapsed time since start()
   Real elapsed() {return _now() - _t0;}

};  // Timer

}   // AMOCH

#endif  // AMOCH_TIMER_H
