// ============================================================================
// chaintype.cc -- AMOCH::Chaintype methods
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#include "chaintype.h"
#include "logfile.h"

using std::string;
using std::vector;
using AMOCH::Chaintype;
using AMOCH::Monomer;
using AMOCH::log_message;
using AMOCH::Param;

// ============================================================================
// Constructor
// ============================================================================
Chaintype::Chaintype(const string& name) :
   _name(name),
   _monomers(0),
   _weights(0),
   _term_monomer(-1),
   _select_by_weight(false)
{}

// ============================================================================
// Return the index of monomers with matching name, or -1 if no match
// ============================================================================
static int 
match_monomer_name(const vector<Monomer>& monomers,
                   const string& name)
{
   int nmon = (int)monomers.size();

   for (int i = 0; i < nmon; i++) {
      if (name == monomers[i].name())
         return i;
   }
   return -1;
}

// ============================================================================
// Set properties from input; return true on success or false on error
// ============================================================================
bool 
Chaintype::setup(const Param& pchtype,
                 const vector<Monomer>& monomers)
{
   const Param *p = pchtype.find("monomers");
   const char *ctname = _name.c_str();

   if (!p) {
      log_message("\nNo monomers specified for chain type %s\n", ctname);
      return false;
   }
   log_message("Chain type %s:\n   monomers", ctname);

   const int nmon = p->nvalues();
   int i, j;
   string mname;

   _monomers.reserve(nmon);
   for (i = 0; i < nmon; i++) {
      mname = p->value(i);
      j = match_monomer_name(monomers, mname);
      if (j == -1) {
         log_message("\nUnknown monomer %s in chain_type %s\n", mname.c_str(), 
                      ctname);
         return false;
      }
      _monomers.push_back(j);
      log_message(" %s", mname.c_str());
   }
   log_message("\n");

   if ((p = pchtype.find("monomer_weights"))) {
      Real wsum = 0.0;

      _select_by_weight = true;
      _weights.reserve(nmon);
      log_message("   weights");
      for (i = 0; i < nmon; i++) {
         _weights.push_back(atof(p->value(i).c_str()));
         log_message(" %g", _weights[i]);
         wsum += _weights[i];
      }
      if (fabs(wsum-1.0) > REAL_SMALL) {
         log_message("\nMonomer weights for chain_type %s do not sum to 1\n",
                     ctname);
         return false;
      }
      log_message("\n");
   }
   else
      log_message("   Pattern repeats\n");

   if ((p = pchtype.find("term"))) {
      mname = p->value();
      j = match_monomer_name(monomers, mname);
      if (j == -1) {
         log_message("\nUnknown terminating monomer %s in chain_type %s\n",
                     mname.c_str(), ctname);
         return false;
      }
      _term_monomer = j;
      log_message("   Terminate chain with monomer %s\n", mname.c_str());
   }
   return true;
}
