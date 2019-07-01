// ============================================================================
// slab.cc -- AMOCH::Slab methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include <cstdlib>
#include "slab.h"
#include "logfile.h"

using std::vector;
using AMOCH::Slab;
using AMOCH::Real;
using AMOCH::Triple;
using AMOCH::log_message;

// ============================================================================
// Constructor
// ============================================================================
Slab::Slab() :
   AMOCH::Exclusion(),
   _corners(8),
   _edges(12),
   _normals(6)
{}

// ============================================================================
// Return true if parameters can be extracted from pexc, else return false
// after logging a message
// ============================================================================
bool 
Slab::setup(const AMOCH::Param& pexc)
{
   const Param *p;
   Triple min, max;

   if ((p = pexc.find("min"))) {
      if (p->nvalues() != 3) {
         log_message("Expecting 3 values [x,y,z] for excluded slab min\n");
         return false;
      }
      min.x = atof(p->value(0).c_str());
      min.y = atof(p->value(1).c_str());
      min.z = atof(p->value(2).c_str());
   }
   else {
      log_message("\nMissing minimum point for excluded slab\n");
      return false;
   }
   if ((p = pexc.find("max"))) {
      if (p->nvalues() != 3) {
         log_message("Expecting 3 values [x,y,z] for excluded slab max\n");
         return false;
      }
      max.x = atof(p->value(0).c_str());
      max.y = atof(p->value(1).c_str());
      max.z = atof(p->value(2).c_str());
   }
   else {
      log_message("\nMissing maximum point for excluded slab\n");
      return false;
   }

   log_message("   Excluded slab from (%g, %g, %g) to (%g, %g, %g)\n", min.x,
               min.y, min.z, max.x, max.y, max.z);

   // 0 == min, 7 == max
   //     
   //             6 ---- 7
   //             |      |
   //  2 ---- 3   |      |   z      positive y into screen
   //  |      |   4 ---- 5   ^
   //  |      |              |
   //  0 ---- 1              |-- >x   

   int i;

   for (i = 0; i < 4; i++)
      _corners[i] = min;
   _corners[1].x = max.x;
   _corners[2].z = max.z;
   _corners[3].x = max.x;
   _corners[3].z = max.z;
   for (i = 4; i < 8; i++)
      _corners[i] = max;
   _corners[4].x = min.x;
   _corners[4].z = min.z;
   _corners[5].z = min.z;
   _corners[6].x = min.x;

   //                to          from   
   _edges[ 0] = _corners[0] - _corners[1];
   _edges[ 1] = _corners[0] - _corners[2];
   _edges[ 2] = _corners[1] - _corners[3];
   _edges[ 3] = _corners[2] - _corners[3];
   _edges[ 4] = _corners[4] - _corners[5];
   _edges[ 5] = _corners[4] - _corners[6];
   _edges[ 6] = _corners[5] - _corners[7];
   _edges[ 7] = _corners[6] - _corners[7];
   _edges[ 8] = _corners[0] - _corners[4];
   _edges[ 9] = _corners[1] - _corners[5];
   _edges[10] = _corners[2] - _corners[6];
   _edges[11] = _corners[3] - _corners[7];
   for (i = 0; i < 12; i++)
      _edges[i].normalize();

   _normals[0] = cross(_edges[0], _edges[2]);   // -y
   _normals[1] = cross(_edges[4], _edges[5]);   // +y
   _normals[2] = cross(_edges[1], _edges[8]);   // -x
   _normals[3] = cross(_edges[6], _edges[9]);   // +x
   _normals[4] = cross(_edges[0], _edges[9]);   // -z
   _normals[5] = cross(_edges[3], _edges[10]);  // +z
   for (i = 0; i < 6; i++)
      _normals[i].normalize();
   return true;
}

// ============================================================================
// Return true if pos is inside the excluded volume, else return false
// ============================================================================
bool 
Slab::inside(const Triple& pos) const
{
   return pos.x > _corners[0].x && 
          pos.y > _corners[0].y && 
          pos.z > _corners[0].z &&
          pos.x < _corners[7].x && 
          pos.y < _corners[7].y && 
          pos.z < _corners[7].z;
}

// ============================================================================
// Set the components of comp to the shortest vector from pos to a point 
// on the indicated edge (_edges[edge]) with vertex at the indicated corner
// ============================================================================
void 
Slab::_comp_to_edge(const Triple& pos,
                    int edge,
                    int corner,
                    Triple& comp) const
{
   Triple q = _edges[edge];
   Triple v = pos - _corners[corner];

   q *= dot(v, q);
   comp = q - v;
}

// ============================================================================
// Set the components of comp to the shortest vector from pos to a point 
// in the indicated face, with normal _normals[face] and _corners[corner]
// in the face plane
// ============================================================================
void 
Slab::_comp_to_face(const Triple& pos,
                    int face,
                    int corner,
                    Triple& comp) const
{
   Triple v = pos - _corners[corner];

   comp = _normals[face] * -dot(v, _normals[face]);
}

// ============================================================================
// Return the minimum square distance from pos to the excluded volume
// ============================================================================
Real 
Slab::min_sqdist(const Triple& pos) const
{
   Triple comp;
   int state_x = 0;
   int state_y = 0;
   int state_z = 0;

   if (pos.x < _corners[0].x)
      state_x = -1;
   if (pos.x > _corners[7].x)
      state_x =  1;
   if (pos.y < _corners[0].y)
      state_y = -1;
   if (pos.y > _corners[7].y)
      state_y =  1;
   if (pos.z < _corners[0].z)
      state_z = -1;
   if (pos.z > _corners[7].z)
      state_z =  1;
   switch (state_x + state_y*10 + state_z*100) {
      case -111:
         comp = _corners[0] - pos;
         break;
      case -110:
         _comp_to_edge(pos, 0, 1, comp);
         break;
      case -109:
         comp = _corners[1] - pos;
         break;
      case -101:
         _comp_to_edge(pos, 8, 4, comp);
         break;
      case -100:
         _comp_to_face(pos, 4, 0, comp);
         break;
      case -99:
         _comp_to_edge(pos, 9, 5, comp);
         break;
      case -91:
         comp = _corners[4] - pos;
         break;
      case -90:
         _comp_to_edge(pos, 4, 5, comp);
         break;
      case -89:
         comp = _corners[5] - pos;
         break;
      case -11:
         _comp_to_edge(pos, 1, 2, comp);
         break;
      case -10:
         _comp_to_face(pos, 0, 0, comp);
         break;
      case -9:
         _comp_to_edge(pos, 2, 3, comp);
         break;
      case -1:
         _comp_to_face(pos, 2, 0, comp);
         break;
      case 1:
         _comp_to_face(pos, 3, 7, comp);
         break;
      case 9:
         _comp_to_edge(pos, 5, 6, comp);
         break;
      case 10:
         _comp_to_face(pos, 1, 7, comp);
         break;
      case 11:
         _comp_to_edge(pos, 6, 7, comp);
         break;
      case 89:
         comp = _corners[2] - pos;
         break;
      case 90:
         _comp_to_edge(pos, 3, 3, comp);
         break;
      case 91:
         comp = _corners[3] - pos;
         break;
      case 99:
         _comp_to_edge(pos, 10, 6, comp);
         break;
      case 100:
         _comp_to_face(pos, 5, 7, comp);
         break;
      case 101:
         _comp_to_edge(pos, 11, 7, comp);
         break;
      case 109:
         comp = _corners[6] - pos;
         break;
      case 110:
         _comp_to_edge(pos, 7, 7, comp);
         break;
      case 111:
         comp = _corners[7] - pos;
         break;
      case 0: // inside the slab...
         comp.zero();
         break;
      default: // ???
         comp.zero();
         log_message("ERROR: Position is in some other dimension...\n");
         break;
   }
   return comp.sqlength();
}

// ============================================================================
// Return the excluded volume in cubic Angstroms
// ============================================================================
Real 
Slab::volume() const
{
   Triple diag = _corners[7] - _corners[0];

   return diag.x * diag.y * diag.z;
}

