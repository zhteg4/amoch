// ============================================================================
// system.h -- AMOCH::System class
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#ifndef AMOCH_SYSTEM_H
#define AMOCH_SYSTEM_H

#include "cell.h"
#include "chain.h"

namespace AMOCH {

//
// Spatial domain of simulation, divided into Cells
//
class System {
   
   private:

   Triple _min;    // minimum point
   Triple _max;    // maximum point
   Triple _size;   // dimensions
   Triple _half;   // half dimensions
   Triple _cell;   // size of a sub-volume ("cell")
   int _nx;        // number of cells in X
   int _ny;        // number of cells in Y
   int _nz;        // number of cells in Z
   int _nxy;       // nx * ny
   int _nxyz;      // nx * ny * nz

   public:

   std::vector<Cell> cells;

   System(const Triple& min,
          const Triple& max,
          Real cell_size,
          int max_particles,   // in simulation
          int max_bonds);      // per particle

   // Synthesized copy constructor, assignment, destructor

   const Triple& min() const {return _min;}
   const Triple& max() const {return _max;}
   const Triple& size() const {return _size;}
   int ncells() const {return _nxyz;}

   // Return the index of the cell that contains the position p
   int cell_index(const Triple& p) const;
          
   // Fill nbr_cells with indices of Cells that surround Cell n
   void find_neighbor_cells(int n,
                            std::vector<int>& nbr_cells) const;

   // Wrap the position p to a point inside
   void wrap(Triple& p) const;
     
   // Return the minimum square distance from posi to posj
   Real min_sqdist(const Triple& posi,
                   const Triple& posj) const;

   // Return true if the particles a and b in the Chain c are separated by 
   // nbonds or fewer bonds, else return false
   bool bond_separation(const Chain& c, 
                        int a,
                        int b,
                        int nbonds) const;
};  // System

}   // AMOCH

#endif  // AMOCH_SYSTEM_H
