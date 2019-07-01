// ============================================================================
// types.h -- AMOCH::PartType, AMOCH::BondType classes
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#ifndef AMOCH_TYPES_H
#define AMOCH_TYPES_H

#include <vector>
#include <string>

namespace AMOCH {

// ----------------------------------------------------------------------------
//
// Particle type info
//
class PartType {

   private:

   int _element;
   std::string _symbol;
   int _nbonds;
   bool _hbond;
   bool _resonant;

   public:

   PartType(int element = 0,
            const std::string& symbol = "",
            int nbonds = -1,
            bool hbond = false,
            bool resonant = false);

   // Synthesized copy constructor, assignment, destructor

   // Readers
   int element() const {return _element;}
   const std::string& symbol() const {return _symbol;}
   int nbonds() const {return _nbonds;}
   bool hbond() const {return _hbond;}
   bool resonant() const {return _resonant;}

};  // PartType

// Return the index of the matching particle type, or -1 if no match
int match_part_type(const std::vector<PartType>& part_types,
                    int element,
                    int nbonds = -1,
                    bool hbond = false,
                    bool resonant = false);

// ----------------------------------------------------------------------------
//
// Bond type info
//
class BondType {

   private:

   int _ptype1;  // part_types index
   int _ptype2;  // part_types index
   int _order;

   public:

   BondType(int ptype1 = 0,
            int ptype2 = 0,
            int order = -1);

   // Synthesized copy constructor, assignment, destructor

   // Readers
   int ptype1() const {return _ptype1;}
   int ptype2() const {return _ptype2;}
   int order() const {return _order;}

};  // BondType

// Return the index of the matching bond type, or -1 if no match
int match_bond_type(const std::vector<BondType>& bond_types,
                    int ptype1,
                    int ptype2,
                    int order = -1);

}   // AMOCH

#endif  // AMOCH_TYPES_H
