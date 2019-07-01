// ============================================================================
// fene.cc -- AMOCH::FENE methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include <cstdlib>
#include <cmath>
#include "fene.h"
#include "logfile.h"

using std::vector;
using AMOCH::FENE;
using AMOCH::Real;
using AMOCH::Param;
using AMOCH::log_message;

// ============================================================================
// Constructor
// ============================================================================
FENE::FENE(int ntypes,
           Real rcut) :
   AMOCH::Energy(ntypes, 4, rcut),
   _params(ntypes, vector<vector<Real> >(ntypes, vector<Real>(10, 0.0)))
{}

// ============================================================================
// Set properties from input; return true on success or false on error
// ============================================================================
bool 
FENE::setup(const Param& pen,
            const vector<AMOCH::PartType>& part_types)
{
   const Param *p, *q;
   Real eps, sig, k, R0;
   const char *sym;
   unsigned i;

   p = pen.find("params");
   while ((p)) {
      sym = p->value().c_str();
      if ((q = p->find("epsilon")))
         eps = atof(q->value().c_str());
      else {
         log_message("\nFENE epsilon not specified for %s\n", sym);
         return false;
      }
      if ((q = p->find("sigma")))
         sig = atof(q->value().c_str());
      else {
         log_message("\nFENE sigma not specified for %s\n", sym);
         return false;
      }
      if ((q = p->find("k")))
         k = atof(q->value().c_str());
      else {
         log_message("\nFENE k not specified for %s\n", sym);
         return false;
      }
      if ((q = p->find("R0")))
         R0 = atof(q->value().c_str());
      else {
         log_message("\nFENE R0 not specified for %s\n", sym);
         return false;
      }
      for (i = 0; i < part_types.size(); i++) {
         if (p->value() == part_types[i].symbol()) {
            _params_in[i][0] = eps;
            _params_in[i][1] = sig;
            _params_in[i][2] = k;
            _params_in[i][3] = R0;
         }
      }
      log_message("   %s: epsilon = %g, sigma = %g, k = %g, R0 = %g\n", sym, 
                  eps, sig, k, R0);
      p = p->find_next("params");
   }
   return true;
}

//_params_in[type][0] = epsilon
//_params_in[type][1] = sigma
//_params_in[type][2] = k
//_params_in[type][3] = R0

// _params[type1][type2][0] = 48 * epsilon * sigma^12
// _params[type1][type2][1] = 24 * epsilon * sigma^6
// _params[type1][type2][2] =  4 * epsilon * sigma^12
// _params[type1][type2][3] =  4 * epsilon * sigma^6
// _params[type1][type2][4] = epsilon
// _params[type1][type2][5] = (r_c)^2 == (2^(1/6)*sigma)^2
// _params[type1][type2][6] = K
// _params[type1][type2][7] = RO^2
// _params[type1][type2][8] = -0.5 * K * R0^2
// _params[type1][type2][9] = 1 / R0^2

// ============================================================================
// Scale interaction parameters
// ============================================================================
void 
FENE::scale(Real f)
{
   int i, j, k;
   Real e, s, s2, s6;
   Real es6, es12;
   Real K, R0, R0sq;

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
         _params[i][j][4] = e;
         _params[i][j][5] = 1.25992104989488*s2;

         K = sqrt(_params_in[i][2]*_params_in[j][2])*f;
         R0 = sqrt(_params_in[i][3]*_params_in[j][3])*f;
         R0sq = R0*R0;
         _params[i][j][6] = K;
         _params[i][j][7] = R0sq;
         _params[i][j][8] = -0.5*K*R0sq;
         _params[i][j][9] = 1.0/R0sq;
      }
   }
   for (i = 0; i < _ntypes; i++) {
      for (j = i+1; j < _ntypes; j++) {
         for (k = 0; k < 10; k++)
            _params[j][i][k] = _params[i][j][k];
      }
   }
}

// ============================================================================
// Calculate the pair energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
FENE::pair_energy(int t1, 
                  int t2,
                  Real rsq) const
{
   Real g = 0.0;

   if (rsq < _params[t1][t2][5]) {
      Real p2 = 1.0/rsq;
      Real p6 = p2*p2*p2;

      g = p6*(p6*_params[t1][t2][2] - _params[t1][t2][3]) + _params[t1][t2][4];
   }
   return g;
}

// ============================================================================
// Calculate the bond energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
FENE::bond_energy(int t1, 
                  int t2,
                  Real rsq) const
{
   Real g = REAL_MAX;

   if (rsq < _params[t1][t2][7])
      g = _params[t1][t2][8] * log(1.0 - rsq*_params[t1][t2][9]);
   return g;
}

// ============================================================================
// Calculate the hard energy between two particles of types t1, t2, 
// separated by square distance rsq
// ============================================================================
Real 
FENE::hard_energy(int t1, 
                  int t2,
                  Real rsq) const
{
   Real p2 = 1.0/rsq;
   Real p6 = p2*p2*p2;
   
   return p6 * p6 * _params[t1][t2][2];
}

