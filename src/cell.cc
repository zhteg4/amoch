// ============================================================================
// cell.cc -- AMOCH::Cell methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include <cstdio>  // FIXME
#include "cell.h"

using std::vector;
using AMOCH::Cell;
using AMOCH::Ivec;
using AMOCH::Triple;

int Cell::nresize = 0;

// ============================================================================
// Constructor
// ============================================================================
Cell::Cell(int max_particles,
           int max_bonds) :
   _max_particles(max_particles),
   _nparticles(0),
   _max_bonds(max_bonds),
   _chain(max_particles, -1),
   _index(max_particles, -1),
   _ptype(max_particles, -1),
   _nbonds(max_particles, 0),
   _bonds(max_particles, vector<Ivec>(max_bonds, Ivec(2))),
   pos(max_particles, Triple()),
   tpos(max_particles, Triple())
{}

// ============================================================================
// Add a new particle; return the slot into which the particle was added
// ============================================================================
int 
Cell::add_particle(int chain,
                   int index,
                   int ptype,
                   const Triple& p,
                   const Triple& tp)
{
   int n = _nparticles;

   if (n == _max_particles) {
      _max_particles = (int)(_max_particles * RESIZE_SCALE);
      if (n == _max_particles)
         _max_particles++;
      _chain.reserve(_max_particles);
      _index.reserve(_max_particles);
      _ptype.reserve(_max_particles);
      _nbonds.reserve(_max_particles);
      _bonds.reserve(_max_particles);
      pos.reserve(_max_particles);
      tpos.reserve(_max_particles);
      for (int i = n; i < _max_particles; i++) {
         _chain.push_back(-1);
         _index.push_back(-1);
         _ptype.push_back(-1);
         _nbonds.push_back(0);
         _bonds.push_back(vector<Ivec>(_max_bonds, Ivec(2)));
         pos.push_back(Triple());
         tpos.push_back(Triple());
      }
      Cell::nresize++;
   }
   _chain[n] = chain;
   _index[n] = index;
   _ptype[n] = ptype;
   pos[n] = p;
   tpos[n] = tp;
   ++_nparticles;
   return n;
}

// ============================================================================
// Add a bond (to particle bindex) to the particle in the indicated slot.
// ============================================================================
void 
Cell::add_bond(int slot,
               int bindex,
               int btype)
{
   int nb = _nbonds[slot];

   if (nb == (int)_bonds[slot].size()) {
      int m = (int)(nb * RESIZE_SCALE);

      if (m == nb)
         m++;
      _bonds[slot].reserve(m);
      for (int i = nb; i < m; i++)
         _bonds[slot].push_back(Ivec(2));
      Cell::nresize++;
   }
   _bonds[slot][nb][0] = bindex;
   _bonds[slot][nb][1] = btype;
   _nbonds[slot]++;
}

// ============================================================================
// Return true if the particle in the indicated slot is bonded to the 
// particle bindex, else return false.
// ============================================================================
bool 
Cell::bonded(int slot,
             int bindex) const
{
   for (int i = 0; i < _nbonds[slot]; i++) {
      if (bindex == _bonds[slot][i][0])
         return true;
   }
   return false;
}

