// ============================================================================
// slab.h -- AMOCH::Slab class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_SLAB_H
#define AMOCH_SLAB_H

#include "exclusion.h"

namespace AMOCH {

//
// Rectangular excluded volume
//
class Slab: public Exclusion {

   private:

   std::vector<Triple> _corners;
   std::vector<Triple> _edges;
   std::vector<Triple> _normals;

   // Set the components of comp to the shortest vector from pos to a point 
   // on the indicated edge (_edges[edge]) with vertex at the indicated corner
   void _comp_to_edge(const Triple& pos,
                      int edge,
                      int corner,
                      Triple& comp) const;

   // Set the components of comp to the shortest vector from pos to a point 
   // in the indicated face, with normal _normals[face] and _corners[corner]
   // in the face plane
   void _comp_to_face(const Triple& pos,
                      int face,
                      int corner,
                      Triple& comp) const;

   public:

   Slab();

   // Synthesized copy constructor, assignment, destructor

   // Return true if parameters can be extracted from pexc, else return false
   // after logging a message
   bool setup(const Param& pexc);

   // Return true if pos is inside the excluded volume, else return false
   bool inside(const Triple& pos) const;

   // Return the minimum square distance from pos to the excluded volume
   Real min_sqdist(const Triple& pos) const;

   // Return the excluded volume in cubic Angstroms
   Real volume() const;

};  // Slab

}   // AMOCH

#endif  // AMOCH_SLAB_H
