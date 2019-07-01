// ============================================================================
// types.cc -- AMOCH::Types methods
// ----------------------------------------------------------------------------
// Copyright 
// ----------------------------------------------------------------------------
// LICENSE 
// ============================================================================

#include "types.h"
#include "logfile.h"

using std::vector;
using AMOCH::PartType;
using AMOCH::BondType;

// ============================================================================
// PartType constructor
// ============================================================================
PartType::PartType(int element,
                   const std::string& symbol,
                   int nbonds,
                   bool hbond,
                   bool resonant) : 
   _element(element),
   _symbol(symbol),
   _nbonds(nbonds),
   _hbond(hbond),
   _resonant(resonant)
{}

// ============================================================================
// Return the index of the matching particle type, or -1 if no match
// ============================================================================
int 
AMOCH::match_part_type(const vector<PartType>& part_types,
                       int element,
                       int nbonds,
                       bool hbond,
                       bool resonant)
{
   int i;

   for (i = 0; i < (int)part_types.size(); i++) {
      if ((element == part_types[i].element()) && 
          (nbonds == part_types[i].nbonds()) &&
          (hbond == part_types[i].hbond()) && 
          (resonant == part_types[i].resonant()))
         return i;
   }
   return -1;
}

// ============================================================================
// BondType constructor
// ============================================================================
BondType::BondType(int ptype1,
                   int ptype2,
                   int order) : 
   _ptype1(ptype1),
   _ptype2(ptype2),
   _order(order)
{}

// ============================================================================
// Return the index of the matching bond type, or -1 if no match
// ============================================================================
int 
AMOCH::match_bond_type(const vector<BondType>& bond_types,
                       int ptype1,
                       int ptype2,
                       int order)
{
   int i;

   for (i = 0; i < (int)bond_types.size(); i++) {
      if ((((ptype1 == bond_types[i].ptype1()) &&
            (ptype2 == bond_types[i].ptype2())) ||
           ((ptype2 == bond_types[i].ptype1()) &&
            (ptype1 == bond_types[i].ptype2()))) && 
          (order == bond_types[i].order()))
         return i;
   }
   return -1;
}

