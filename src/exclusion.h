// ============================================================================
// exclusion.h -- AMOCH::Exclusion class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_EXCLUSION_H
#define AMOCH_EXCLUSION_H

#include "param.h"
#include "triple.h"

namespace AMOCH {

//
// Generic excluded volume
//
// Abstract base class for Sphere, Cylinder, Slab
//
class Exclusion {

   public:

   Exclusion() {}

   // Synthesized copy constructor, assignment

   virtual ~Exclusion() {}

   // Return true if parameters can be extracted from pexc, else return false
   // after logging a message
   virtual bool setup(const Param& pexc) = 0;

   // Return true if pos is inside the excluded volume, else return false
   // XXX this should really be "contains()" or something like that ...
   virtual bool inside(const Triple& pos) const = 0;

   // Return the minimum square distance from pos to the excluded volume
   virtual Real min_sqdist(const Triple& pos) const = 0;

   // Return the excluded volume in cubic Angstroms
   virtual Real volume() const = 0;

};  // Exclusion

}   // AMOCH

#endif  // AMOCH_EXCLUSION_H
