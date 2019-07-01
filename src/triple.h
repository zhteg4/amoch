// ============================================================================
// triple.h -- AMOCH::Triple class
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_TRIPLE_H
#define AMOCH_TRIPLE_H

#include <cmath>
#include "config.h"

namespace AMOCH {

//
// Cartesian triple (X,Y,Z) 
//
class Triple {

   public:

   Real x;
   Real y;
   Real z;

   Triple(Real xx = 0.0,
          Real yy = 0.0,
          Real zz = 0.0);

   // Synthesized copy constructor, assignment, destructor

   // Add-assign
   Triple& operator+=(const Triple& rhs);

   // Subtract-assign
   Triple& operator-=(const Triple& rhs);

   // Scale-assign
   Triple& operator*=(Real f);

   // Return the square (Euclidean) length 
   Real sqlength() const {return x*x + y*y + z*z;}

   // Return the (Euclidean) length 
   Real length() const {return sqrt(x*x + y*y + z*z);}

   // Scale components to unit length
   void normalize() {*this *= 1.0/length();}

   // Set all components to zero
   void zero() {x = y = z = 0.0;}

};  // Triple

}   // AMOCH

// Return lhs + rhs
AMOCH::Triple operator+(const AMOCH::Triple& lhs,
                        const AMOCH::Triple& rhs);

// Return lhs - rhs
AMOCH::Triple operator-(const AMOCH::Triple& lhs,
                        const AMOCH::Triple& rhs);

// Return lhs * f 
AMOCH::Triple operator*(const AMOCH::Triple& lhs,
                        AMOCH::Real f);

// Return the dot (inner) product lhs * rhs
AMOCH::Real dot(const AMOCH::Triple& lhs,
                const AMOCH::Triple& rhs);

// Return the cross (vector) product lhs x rhs
AMOCH::Triple cross(const AMOCH::Triple& lhs,
                    const AMOCH::Triple& rhs);

#endif  // AMOCH_TRIPLE_H
