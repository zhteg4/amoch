// ============================================================================
// system.cc -- AMOCH::System methods
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#include <cmath>
#include "system.h"

// Periodic edges
#define XMIN  0
#define XMAX  1
#define YMIN  2
#define YMAX  3
#define ZMIN  4
#define ZMAX  5

using std::vector;
using AMOCH::System;
using AMOCH::Triple;
using AMOCH::Real;

// ============================================================================
// Constructor
// ============================================================================
System::System(const Triple& min,
               const Triple& max,
               Real cell_size,
               int max_particles,
               int max_bonds) :
   _min(min), 
   _max(max),
   _size(max-min),
   _half(_size*0.5), 
   _cell(Triple(cell_size, cell_size, cell_size)),
   _nx((int)ceil(_size.x/_cell.x)),
   _ny((int)ceil(_size.y/_cell.y)),
   _nz((int)ceil(_size.z/_cell.z)),
   _nxy(_nx*_ny),
   _nxyz(_nxy*_nz),
#define NEXT2(v) (int)pow(2.0, ceil(log2(v))) 
   cells(_nxyz, AMOCH::Cell(NEXT2((Real)max_particles/(Real)_nxyz), max_bonds))
{}

// ============================================================================
// Return the index of the cell that contains the position p
// ============================================================================
int 
System::cell_index(const Triple& p) const
{
   return (int)((p.x - _min.x) / _cell.x) +
          (int)((p.y - _min.y) / _cell.y) * _nx +
          (int)((p.z - _min.z) / _cell.z) * _nxy;
}
          
// ============================================================================
// Fill nbr_cells with indices of Cells that surround Cell n
// ============================================================================
void 
System::find_neighbor_cells(int n,
                            vector<int>& nbr_cells) const
{
   const int bz = n / _nxy;
   const int by = (n - bz*_nxy) / _nx;
   const int bx = n % _nx;
   const int i_shift[26][3] = {
      { 1,  0,  0}, { 1,  1,  0}, { 0,  1,  0}, {-1,  1,  0}, {-1,  0,  0},
      {-1, -1,  0}, { 0, -1,  0}, { 1, -1,  0}, { 0,  0, -1}, { 1,  0, -1},
      { 1,  1, -1}, { 0,  1, -1}, {-1,  1, -1}, {-1,  0, -1}, {-1, -1, -1},
      { 0, -1, -1}, { 1, -1, -1}, { 0,  0,  1}, { 1,  0,  1}, { 1,  1,  1},
      { 0,  1,  1}, {-1,  1,  1}, {-1,  0,  1}, {-1, -1,  1}, { 0, -1,  1},
      { 1, -1,  1}
   };
   bool edge[6];
   int i;

   edge[XMIN] = (bx == 0);
   edge[XMAX] = (bx == _nx - 1);
   edge[YMIN] = (by == 0);
   edge[YMAX] = (by == _ny - 1);
   edge[ZMIN] = (bz == 0);
   edge[ZMAX] = (bz == _nz - 1);

   // Assume n is not on an edge, initially
   for (i = 0; i < 26; i++)
      nbr_cells[i] = n + i_shift[i][0] + i_shift[i][1]*_nx + i_shift[i][2]*_nxy;

   const int j_xmin[9] = {3, 4, 5, 12, 13, 14, 21, 22, 23};
   const int j_xmax[9] = {0, 1, 7,  9, 10, 16, 18, 19, 25};
   const int j_ymin[9] = {5, 6, 7, 14, 15, 16, 23, 24, 25};
   const int j_ymax[9] = {1, 2, 3, 10, 11, 12, 19, 20, 21};

   // Periodic boundaries
   if (edge[XMIN]) {
      for (i = 0; i < 9; i++)
         nbr_cells[j_xmin[i]] += _nx;
   }
   if (edge[XMAX]) {
      for (i = 0; i < 9; i++)
         nbr_cells[j_xmax[i]] -= _nx;
   }
   if (edge[YMIN]) {
      for (i = 0; i < 9; i++)
         nbr_cells[j_ymin[i]] += _nxy;
   }
   if (edge[YMAX]) {
      for (i = 0; i < 9; i++)
         nbr_cells[j_ymax[i]] -= _nxy;
   }
   if (edge[ZMIN]) {
      for (i = 8; i < 17; i++)
         nbr_cells[i] += _nxyz;
   }
   if (edge[ZMAX]) {
      for (i = 17; i < 26; i++)
         nbr_cells[i] -= _nxyz;
   }

   // Set any duplicate indices to -1
   for (i = 0; i < 26; i++) {
      for (int j = i+1; j < 26; j++) {
         if (nbr_cells[j] == nbr_cells[i]) {
            nbr_cells[j] = -1;
            break;
         }
      }
      if (nbr_cells[i] == n)
         nbr_cells[i] = -1;
   }
}

// ============================================================================
// Wrap the position p to a point inside
// ============================================================================
void 
System::wrap(Triple& p) const
{
   while (p.x >= _max.x)
      p.x -= _size.x;
   while (p.x <  _min.x)
      p.x += _size.x;
   while (p.y >= _max.y)
      p.y -= _size.y;
   while (p.y <  _min.y)
      p.y += _size.y;
   while (p.z >= _max.z)
      p.z -= _size.z;
   while (p.z <  _min.z)
      p.z += _size.z;
}
     
// ============================================================================
// Return the minimum square distance from posi to posj, and store the 
// displacement components in comp
// ============================================================================
Real 
System::min_sqdist(const Triple& posi,
                   const Triple& posj) const
{
   Triple comp = posj - posi;

   if (comp.x > _half.x)
      comp.x -= _size.x;
   if (comp.x < -_half.x)
      comp.x += _size.x;
   if (comp.y > _half.y)
      comp.y -= _size.y;
   if (comp.y < -_half.y)
      comp.y += _size.y;
   if (comp.z > _half.z)
      comp.z -= _size.z;
   if (comp.z < -_half.z)
      comp.z += _size.z;
   return comp.sqlength();
}

// ============================================================================
// Return true if the particles a and b in the Chain c are separated by nbonds
// or fewer bonds, else return false
//
// Breadth-first search out to nbonds edges in the bond graph with source a
// ============================================================================
bool
System::bond_separation(const AMOCH::Chain& c, 
                        int a,
                        int b,
                        int nbonds) const
{
   const int a_cell = c.cell(a);
   const int a_slot = c.slot(a);
   const int nba = cells[a_cell].nbonds(a_slot);
   vector<AMOCH::Ivec> ba(4, AMOCH::Ivec(2));

   if (nbonds == 0)
      return false;
   ba = cells[a_cell].bonds(a_slot);
   for (int i = 0; i < nba; i++) {
      if (ba[i][0] == b)
         return true;
      if (bond_separation(c, ba[i][0], b, nbonds-1))
         return true;
   }
   return false;
}

