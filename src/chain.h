// ============================================================================
// chain.h -- AMOCH::Chain class
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_CHAIN_H
#define AMOCH_CHAIN_H

#include <vector>
#include <cstdio>
#include "triple.h"

namespace AMOCH {

//
// Molecule built from repeating Monomers; note this class stores info about
// the molecule, but not the particles (see Cell for particles)
//
class Chain {

   private:

   int _type;
   int _max_monomers;
   int _max_particles;
   int _nmonomers;                   // current number of monomers
   int _nparticles;                  // current number of particles
   std::vector<int> _monomer_start;  // index of first particle in each monomer
   std::vector<int> _torsion_count;  // torsion angle selections
   std::vector<Real> _length;        // length as a function of monomers
   std::vector<Triple> _com;         // monomer center of mass
   std::vector<Triple> _m0;          // position of first monomer particle
   std::vector<int> _cell;           // cell index of each particle
   std::vector<int> _slot;           // cell slot of each particle

   public:

   int tail_index;                   // index of current tail particle

   Chain(int type = -1,
         int max_monomers = 0,
         int max_particles = 0);

   // Synthesized copy constructor, assignment, destructor

   // Readers
   int type() const {return _type;}
   int nmonomers() const {return _nmonomers;}
   int max_monomers() const {return _max_monomers;}
   int nparticles() const {return _nparticles;}
   int torsion_count(int n) const {return _torsion_count[n];}
   Real length(int n) const {return _length[n];}
   const Triple& com(int n) const {return _com[n];}
   int cell(int n) const {return _cell[n];}
   int slot(int n) const {return _slot[n];}

   // Add cell, slot info for particle n
   void add_particle_cell(int n,
                          int cell,
                          int slot);

   // Add a monomer consisting of nparticles to the Chain
   void add_monomer(int nparticles);

   // Record a final torsion angle
   void record_torsion(Real angle);

   // Record a final particle position, after all building is complete
   void record_particle(int index,
                        const Triple& pos);

   // Calculate final chain statistics
   void finalize(const Triple& tail_position);

   // Write internal state to the output FILE f
   void write(FILE *f) const;

   // Read internal state from the input FILE f
   void read(FILE *f);

};  // Chain

}   // AMOCH

#endif  // AMOCH_CHAIN_H
