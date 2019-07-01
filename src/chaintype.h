// ============================================================================
// chaintype.h -- AMOCH::Chaintype class
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_CHAINTYPE_H
#define AMOCH_CHAINTYPE_H

#include <string>
#include <vector>
#include "monomer.h"
#include "param.h"
#include "rng.h"

namespace AMOCH {

//
// Monomer arrangement used to build a Chain
//
class Chaintype {

   private:

   std::string _name;
   std::vector<int> _monomers;  // indices
   std::vector<Real> _weights;  // monomer probabilities
   int _term_monomer;           // index
   bool _select_by_weight;

   public:

   explicit Chaintype(const std::string& name = "");

   // Synthesized copy constructor, assignment, destructor

   // Set properties from input; return true on success or false on error
   bool setup(const Param& pchtype,
              const std::vector<Monomer>& monomers);

   // Property readers
   const std::string& name() const {return _name;}
   int nmonomers() const {return _monomers.size();}
   int monomer_index(int n) const {return _monomers[n];}
   Real monomer_weight(int n) const {return _weights[n];}
   bool select_by_weight() const {return _select_by_weight;}
   int monomer_index_by_weight(RNG& rng) const {
      return rng.select_by_weight(_weights);
   }
   int term_monomer_index() const {return _term_monomer;}

};  // Chaintype

}   // AMOCH

#endif  // AMOCH_CHAINTYPE_H
