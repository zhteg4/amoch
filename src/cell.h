// ============================================================================
// cell.h -- AMOCH::Cell class
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_CELL_H
#define AMOCH_CELL_H

#include <vector>
#include "ivec.h"
#include "triple.h"

namespace AMOCH {

//
// Small volume with particles
//
class Cell {

   private:

   int _max_particles;
   int _nparticles;
   int _max_bonds;
   // Per-particle data
   std::vector<int> _chain;
   std::vector<int> _index;   // within Chain
   std::vector<int> _ptype;   // PartType index
   std::vector<int> _nbonds;
   std::vector<std::vector<Ivec> > _bonds;  // [particle][bond][index,btype]

   public:

   std::vector<Triple> pos;   // periodic (wrapped) position
   std::vector<Triple> tpos;  // true (unwrapped) position
   static int nresize;        // number of resize operations

   Cell(int max_particles,
        int max_bonds);       // per particle

   // Synthesized copy constructor, assignment, destructor

   // Readers for slot n
   int nparticles() const {return _nparticles;}
   int chain(int n) const {return _chain[n];}
   int index(int n) const {return _index[n];}
   int ptype(int n)  const {return _ptype[n];}
   int nbonds(int n) const {return _nbonds[n];}
   const std::vector<Ivec>& bonds(int n) const {return _bonds[n];}

   // Add a new particle; return the slot into which the particle was added
   int add_particle(int chain,
                    int index,
                    int ptype,
                    const Triple& p,
                    const Triple& tp);

   // Add a bond (to particle bindex) to the particle in the indicated slot.
   void add_bond(int slot,
                 int bindex,
                 int btype);

   // Return true if the particle in the indicated slot is bonded to the 
   // particle bindex, else return false.
   bool bonded(int slot,
               int bindex) const;

};  // Cell

}   // AMOCH

#endif  // AMOCH_CELL_H
