// ============================================================================
// monomer.h -- AMOCH::Monomer class
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_MONOMER_H
#define AMOCH_MONOMER_H

#include <cstdio>
#include <openbabel/mol.h>
#include "triple.h"
#include "types.h"
#include "ivec.h"
#include "param.h"

namespace AMOCH {

// ----------------------------------------------------------------------------
//
// Internal coordinates for one atom
//
struct ZEntry {
   int type;            // types index
   int bond_atom;
   Real bond_length;    // Angstrom
   int angle_atom;
   Real bond_angle;     // degrees
   int torsion_atom;
   Real torsion_angle;  // degrees

   ZEntry(int typ = -1,
          int bat = -1,
          Real blen = 0.0,
          int aat = -1,
          Real bang = 0.0,
          int tat = -1,
          Real tang = 0.0);

   // Synthesized copy constructor, assignment, destructor

};  // ZEntry

// ----------------------------------------------------------------------------
//
// Repeating molecule in a chain
//
class Monomer {

   private:

   std::string _name;
   OpenBabel::OBMol _mol;
   int _natoms;
   int _nbb;                      // number of backbone atoms
   int _nbonds;
   Triple _anchor;                // position used to fix atom 2
   Real _mass;                    // molecule mass (amu)
   Real _bb_mass;                 // backbone mass (amu)
   Real _length;                  // max distance between particles (Ang)
   std::vector<Ivec> _bonds;      // all bonds in molecule
   std::vector<ZEntry> _zmatrix;  // internal coordinates

   // Add non-backbone atoms to id and reorder molecule
   void _reorder(std::vector<int> id);

   public:

   std::vector<Triple> position;  // atomic positions 

   explicit Monomer(const std::string& mname = "");

   // Synthesized copy constructor, assignment, destructor

   // Read molecule from file.  Return true on success, false on error with 
   // message written to log file.
   bool read(const char *path);

   // Identify backbone atoms from head to tail;  reorder molecule so backbone 
   // atoms are first, followed by all side group atoms.
   void identify_backbone(int head,
                          int tail);

   // Specify backbone atoms from head to tail;  reorder molecule so backbone 
   // atoms are first, followed by all side group atoms.
   void specify_backbone(const std::vector<int>& backbone);

   // Set iternal coordinates, along with positions, mass, and length; 
   // add new particle types to part_types.  Return true on success, or false
   // on error after logging a message.
   bool setup(const Param& pmon,
              std::vector<PartType>& part_types);

   // Identify bonds; if internal coordinates are set, report extra bonds not 
   // in _zmatrix.  Add new bond types to bond_types
   void find_bonds(std::vector<BondType>& bond_types);

   // Adjust internal coordinates to ensure torsion rotations preserve 
   // structure
   void adjust_internal_coords(const std::vector<PartType>& part_types);

   // Property readers
   const std::string& name() const {return _name;}
   int natoms() const {return _natoms;}
   int nbackbone() const {return _nbb;}
   int nbonds() const {return _nbonds;}
   Real mass() const {return _mass;}
   Real bbmass() const {return _bb_mass;}
   Real length() const {return _length;}
   const Ivec& bond(int n) const {return _bonds[n];}

   // Atom n property readers
   int type(int n) const {return _zmatrix[n].type;}
   Real bond_length(int n) const {return _zmatrix[n].bond_length;}
   Real torsion_angle(int n) const {return _zmatrix[n].torsion_angle;}

   // Set the torsion associated with atom n to angle (degrees)
   void set_torsion(int n,
                    Real angle) 
   {
      _zmatrix[n].torsion_angle = angle;
   }

   // Set anchor position
   void set_anchor(const Triple& pos) {_anchor = pos;}

   // Update atomic positions from internal coordinates; if tail_length > 0.0,
   // use it as the length of the tail atom bond.  
   //
   // Requires position[0], position[1], and _anchor
   void update_positions(Real tail_length = 0.0);

   // Return the number of bonds between atoms m and n (m < n); returns _nbonds
   // if no path is found (error)
   int bond_separation(int m,
                       int n) const;

   // Return the maximum number of bonds on any atom in the Monomer m
   int max_bonds_per_atom() const;

   // Write the zmatrix to the output FILE f
   void write_zmatrix(FILE *f,
                      const std::vector<PartType>& part_types) const;

   // Write atomic postitions to an output file
   void write(const char *path);

};  // Monomer

}   // AMOCH

#endif  // AMOCH_MONOMER_H
