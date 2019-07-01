// ============================================================================
// lj.cc -- AMOCH::LJ methods
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#include <cstdlib>
#include <cmath>
#include "lj.h"
#include "logfile.h"

using std::vector;
using AMOCH::LJ;
using AMOCH::log_message;
using AMOCH::Real;
using AMOCH::Param;

// ============================================================================
// Constructor
// ============================================================================
LJ::LJ(int ntypes,
       Real rcut) : 
   AMOCH::Energy(ntypes, 2, rcut),
   _params(ntypes, vector<vector<Real> >(ntypes, vector<Real>(4, 0.0)))
{}

// ============================================================================
// Set properties from input; return true on success or false on error
// ============================================================================
bool 
LJ::setup(const Param& pen,
          const vector<AMOCH::PartType>& part_types)
{
   const Param *p, *q;
   Real eps, sig;
   const char *sym;
   unsigned i;

   p = pen.find("params");
   while ((p)) {
      sym = p->value().c_str();
      if ((q = p->find("epsilon")))
         eps = atof(q->value().c_str());
      else {
         log_message("\nLJ epsilon not specified for %s\n", sym);
         return false;
      }
      if ((q = p->find("sigma")))
         sig = atof(q->value().c_str());
      else {
         log_message("\nLJ sigma not specified for %s\n", sym);
         return false;
      }
      for (i = 0; i < part_types.size(); i++) {
         if (p->value() == part_types[i].symbol()) {
            _params_in[i][0] = eps;
            _params_in[i][1] = sig;
         }
      }
      log_message("   %s: epsilon = %g, sigma = %g\n", sym, eps, sig);
      p = p->find_next("params");
   }
   return true;
}

//_params_in[type][0] = epsilon
//_params_in[type][1] = sigma

// _params[type1][type2][0] = 48 * epsilon * sigma^12
// _params[type1][type2][1] = 24 * epsilon * sigma^6
// _params[type1][type2][2] =  4 * epsilon * sigma^12
// _params[type1][type2][3] =  4 * epsilon * sigma^6

// ============================================================================
// Scale interaction parameters
// ============================================================================
void 
LJ::scale(Real f)
{
   int i, j, k;
   Real e, s, s2, s6;
   Real es6, es12;

   for (i = 0; i < _ntypes; i++) {
      for (j = i; j < _ntypes; j++) {
         e = sqrt(_params_in[i][0]*_params_in[j][0])*f;
         s = sqrt(_params_in[i][1]*_params_in[j][1])*f;
         s2 = s*s;
         s6 = s2*s2*s2;
         es6 = e*s6;
         es12 = es6*s6;
         _params[i][j][0] = 48.0*es12;
         _params[i][j][1] = 24.0*es6;
         _params[i][j][2] = 4.0*es12;
         _params[i][j][3] = 4.0*es6;
      }
   }
   for (i = 0; i < _ntypes; i++) {
      for (j = i+1; j < _ntypes; j++) {
         for (k = 0; k < 4; k++) {
            _params[j][i][k] = _params[i][j][k];
         }
      }
   }
}

// ============================================================================
// Calculate the pair energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
LJ::pair_energy(int t1, 
                int t2,
                Real rsq) const
{
   Real p2 = 1.0/rsq;
   Real p6 = p2*p2*p2;

   return p6 * (p6 * _params[t1][t2][2] - _params[t1][t2][3]);
}

// ============================================================================
// Calculate the bond energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
LJ::bond_energy(int t1, 
                int t2,
                Real rsq) const
{
   return 0.0 * rsq * _params[t1][t2][2];
}

// ============================================================================
// Calculate the hard energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
LJ::hard_energy(int t1, 
                int t2,
                Real rsq) const
{
   Real p2 = 1.0/rsq;
   Real p6 = p2*p2*p2;

   return p6 * p6 * _params[t1][t2][2];
}

