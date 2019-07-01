// ============================================================================
// cbmc.cc -- AMOCH::CBMC methods
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#include <cstdlib>
#include <cmath>
#include "cbmc.h"
#include "logfile.h"
#include "defs.h"  // kB
#include "output.h"

using std::vector;
using AMOCH::CBMC;
using AMOCH::Triple;
using AMOCH::Monomer;
using AMOCH::Param;
using AMOCH::log_message;
using AMOCH::BondType;
using AMOCH::PartType;
using AMOCH::match_part_type;
using AMOCH::match_bond_type;
using AMOCH::System;
using AMOCH::Chain;
using AMOCH::Exclusion;
using AMOCH::Energy;
using AMOCH::RNG;
using AMOCH::Real;
using AMOCH::Ivec;

enum TorsionOption {
   TORSION_SP3,
   TORSION_RANDOM,
   TORSION_CUSTOM
};

// ============================================================================
// Constructor
// ============================================================================
CBMC::CBMC() :
   AMOCH::Builder(),
   _nbuild(-1),
   _nconfigs(30),
   _ndelta_steps(0),
   _bonded_cutoff(4),
   _nchains(0),
   _bb_bondtype(-1),
   _bb_only(false),
   _rotate(true),
   _only_between(false),
   _log_monomers(false),
   _log_energies(false),
   _log_status(false),
   _inv_kT(0.0),
   _delta_degs(3.0),    // degrees
   _bb_bondlen(1.53),   // Angstroms
   _torsion_angle(0, 0.0),
   _torsion_energy(0, 0.0),
   _torsion_prob(0, 0.0),
   _pattern_index(0),
   _chain_tail_anchor(0),
   _chain_tail_pos(0),
   _chain_tail_vec(0)
{}

// ============================================================================
// Setup Monomer m with parameters in pmon; return true on success or 
// false on error after logging a message
// ============================================================================
bool 
CBMC::setup_monomer(Monomer& m,
                    const Param& pmon)
{
   int head, tail;
   const char *mname = m.name().c_str();
   const Param *p;

   if ((p = pmon.find("backbone"))) {
      vector<int> bb(p->nvalues(), -1);

      for (int i = 0; i < p->nvalues(); i++)
         bb[i] = atoi(p->value(i).c_str());
      m.specify_backbone(bb);
   }
   else  {
      if ((p = pmon.find("head"))) 
         head = atoi(p->value().c_str()) - 1;
      else {
         log_message("\nHead atom not specified for monomer %s\n", mname);
         return false;
      }

      if ((p = pmon.find("tail"))) 
         tail = atoi(p->value().c_str()) - 1;
      else {
         log_message("\nTail atom not specified for monomer %s\n", mname);
         return false;
      }
      m.identify_backbone(head, tail);
   }
   return true;
}
                              
// ============================================================================
// Update types of particles and bonds, if any changes when forming bonds
// between monomers m1 and m2
// ============================================================================
void 
CBMC::update_types(const Monomer& m1,
                   const Monomer& m2,
                   vector<PartType>& part_types,
                   vector<BondType>& bond_types)
{
   // No changes to part_types
   if (0 == part_types.size())
      log_message("\nNo particle types!\n");

   const int order = 1;
   const int ptype1 = m1.type(m1.nbackbone()-2);  // m1 tail neighbor type
   const int ptype2 = m2.type(1);                 // m2 head neighbor type
   int btype = match_bond_type(bond_types, ptype1, ptype2, order);

   if (-1 == btype)
      bond_types.push_back(BondType(ptype1, ptype2, order));
}

// ============================================================================
// Set properties from input; return true on success or false on error after
// logging a message
// ============================================================================
bool 
CBMC::setup(const Param& pbuild,
            const vector<PartType>& part_types,
            vector<BondType>& bond_types)
{
   const Param *p, *q;
   Real T = 300.0;

   log_message("\nBuild chains using the Configurational Bias ");
   log_message("Monte Carlo method\n");

   // Build ? monomers
   if ((p = pbuild.find("build")))  {
      _nbuild = atoi(p->value().c_str());
      log_message("   Add %d monomers to chains in the system\n", _nbuild);
   }
   else
      log_message("   Add all monomers to all chains in the system\n");

   // Rotate to configure?
   if ((p = pbuild.find("rotate_torsions"))) {
      if (p->value() == "all") {
         _rotate = true;
         _only_between = false;
      }
      else if (p->value() == "none")
         _rotate = false;
      else if (p->value() == "between_monomers") {
         _rotate = true;
         _only_between = true;
      }
      else {
         log_message("\nUnknown argument \"%s\" to rotate_torsions\n", 
                     p->value().c_str());
         return false;
      }
   }
   if (_rotate) {
      log_message("   Rotate torsions to configure monomers\n");
      if (_only_between)
         log_message("      Rotate only torsions between monomers\n");
   }

   // Number of MC configurations per torsion
   if ((p = pbuild.find("configs")))
      _nconfigs = atoi(p->value().c_str());
   log_message("   Consider %d configurations for each torsion\n", _nconfigs);

   // 1/(kT)
   if ((p = pbuild.find("temperature")))
      T = atof(p->value().c_str());
   log_message("   Build at %g K\n", T);
   _inv_kT = 1.0 / (kB*T);

   // Torsion angle-energy options
   TorsionOption toropt(TORSION_SP3);

   if ((p = pbuild.find("torsion_probs"))) {
      if (p->value() == "sp3")
         toropt = TORSION_SP3;
      else if (p->value() == "random")
         toropt = TORSION_RANDOM;
      else if (p->value() == "custom")
         toropt = TORSION_CUSTOM;
      else {
         log_message("\nUnknown \"torsions\" value: %s\n", p->value().c_str());
         return false;
      }
   }

   Real ang;
   int i, n;

   // Set _torsion_angle, allowed angles, _torsion_energy, and _torsion_prob
   switch (toropt) {
      case TORSION_SP3:
         log_message("   Set torsions using sp3 angle-energies\n");
         _torsion_angle.reserve(360);
         _torsion_energy.reserve(360);
         _torsion_prob.reserve(360);
         for (i = 0; i < 360; i++) {
            ang = (Real)i;
            _torsion_angle.push_back(ang);
            // See DREIDING paper, Mayo et al., J. Phys. Chem. 94, 8897, 1990,
            // eqs (13) and (14)
            _torsion_energy.push_back(1.0 - cos(DEG2RAD(3.0*(ang-180.0))));
         }
         break;
      case TORSION_RANDOM:
         log_message("   Set torsions to random values with ");
         log_message("uniform probability\n");
         _torsion_angle.reserve(360);
         _torsion_energy.reserve(360);
         _torsion_prob.reserve(360);
         for (i = 0; i < 360; i++) {
            _torsion_angle.push_back((Real)i);
            _torsion_energy.push_back(0.0);
         }
         break;
      case TORSION_CUSTOM:
         log_message("   Set torsions using custom angle-energy values\n");
         if ((p = pbuild.find("torsion_angle"))) {
            _torsion_angle.reserve(p->nvalues());
            for (i = 0; i < p->nvalues(); i++) 
               _torsion_angle.push_back(atof(p->value(i).c_str()));
         }
         else {
            log_message("\nMissing torsion_angle values for custom ");
            log_message("angle-energy values\n");
            return false;
         }
         if ((p = pbuild.find("torsion_energy"))) {
            n = p->nvalues();
            if (n != (int)_torsion_angle.size()) {
               log_message("\ntorsion_angle has %d values, but ", 
                           _torsion_angle.size());
               log_message("torsion_energy has %d values\n", n);
               return false;
            }
            _torsion_energy.reserve(n);
            for (i = 0; i < n; i++)
               _torsion_energy.push_back(atof(p->value(i).c_str()));
            _torsion_prob.reserve(n);
         }
         else {
            log_message("\nMissing torsion_energy values for custom ");
            log_message("angle-energy values\n");
            return false;
         }
         break;
   }

   Real Emin = REAL_MAX;
   Real psum = 0.0;

   // Calculate (normalized) probabilities based on energies 
   n = _torsion_energy.size();
   for (i = 0; i < n; i++) {
      if (_torsion_energy[i] < Emin)
         Emin = _torsion_energy[i];
   }
   for (i = 0; i < n; i++) {
      _torsion_prob[i] = exp((Emin - _torsion_energy[i]) * _inv_kT);
      psum += _torsion_prob[i];
   }
   if (fabs(psum-1.0) > REAL_SMALL) {
      for (i = 0; i < n; i++)
         _torsion_prob[i] *= 1.0/psum;
   }

   /*
   log_message("\n\nAngle - probability\n");
   for (i = 0; i < n; i++)
      log_message("%g  %g\n", _torsion_angle[i], _torsion_prob[i]);
   log_message("\n\n");
   */

   // Extra torsion steps?
   if ((p = pbuild.find("extra_steps")))
      _ndelta_steps = atoi(p->value().c_str());
   if ((p = pbuild.find("extra_step_size")))
      _delta_degs = atof(p->value().c_str());
   log_message("   Attempt %d additional steps of +/- %g degrees\n", 
               _ndelta_steps, _delta_degs);

   // Backbone bond; C sp3 - C sp3 by default
   int element1 = 6;
   int element2 = 6;
   int nbonds1 = 4;
   int nbonds2 = 4;
   int order = 1;
   int ptype1, ptype2;

   if ((p = pbuild.find("backbone_bond"))) {
      _bb_bondlen = atof(p->value().c_str());
      if ((q = p->find("element1")))
         element1 = atoi(q->value().c_str());
      if ((q = p->find("nbonds1")))
         nbonds1 = atoi(q->value().c_str());
      if ((q = p->find("element2")))
         element2 = atoi(q->value().c_str());
      if ((q = p->find("nbonds2")))
         nbonds2 = atoi(q->value().c_str());
      if ((q = p->find("order")))
         order = atoi(q->value().c_str());
   }
   log_message("   Bond length between monomers: %g Angstroms\n", _bb_bondlen);
   log_message("   Bond type between monomers: %d with %d bonds - ", 
               element1, nbonds1);
   log_message("%d with %d bonds, order %d\n", element2, nbonds2, order);
   ptype1 = match_part_type(part_types, element1, nbonds1);
   if (-1 == ptype1) {
      log_message("\nUnknown particle type for backbone bond: %d with %d bonds\n",
                  element1, nbonds1);
      return false;
   }
   ptype2 = match_part_type(part_types, element2, nbonds2);
   if (-1 == ptype2) {
      log_message("\nUnknown particle type for backbone bond: %d with %d bonds\n",
                  element2, nbonds2);
      return false;
   }
   _bb_bondtype = match_bond_type(bond_types, ptype1, ptype2, order);
   if (-1 == _bb_bondtype) {
      _bb_bondtype = bond_types.size();
      log_message("      Adding new bond type: %d\n", _bb_bondtype+1);
      bond_types.push_back(BondType(ptype1, ptype2, order));
   }
   else
      log_message("      Bond type: %d\n", _bb_bondtype+1);

   // Nonbonded cutoff
   if ((p = pbuild.find("bonded_cutoff")))
      _bonded_cutoff = atoi(p->value().c_str());
   log_message("   Nonbonded interactions require at least %d bonds ", 
               _bonded_cutoff);
   log_message("between particles\n");

   // Backbone only?
   if ((p = pbuild.find("only"))) {
      if ("backbone" == p->value()) {
         _bb_only = true;
         log_message("   Build only chain backbones, not side groups\n");
      }
      else {
         log_message("\nUnknown \"only\" argument: %s\n", p->value().c_str());
         return false;
      }
   }

   // Log while building?
   p = pbuild.find("log");
   while ((p)) {
      if (p->value() == "monomers")
         _log_monomers = true;
      else if (p->value() == "energies")
         _log_energies = true;
      else if (p->value() == "status")
         _log_status = true;
      else {
         log_message("\nUnknown \"report\" option: %s\n", p->value().c_str());
         return false;
      }
      p = p->find_next("log");
   }
   return true;
}

// ============================================================================
// Return the sum of nonbonded pair interactions of particles in m with
// surrounding particles in system and exclusions
// ============================================================================
Real 
CBMC::_nonbond_interactions(Monomer& m, 
                            const vector<int>& index_map,
                            const Chain& c,
                            int ic,
                            int itail,
                            const System& sys, 
                            const vector<const Exclusion *>& exclusions,
                            const Energy& energy) const
{
   Real E = 0.0;
   Real rsq;
   int i, ti;
   int j, tj;
   int icell, k;
   int m_end, j_ndx;
   int cell_bonds;
   unsigned iexc;
   Triple r, pos;
   vector<int> nbr_cells(26, -1);

   // The index of the first backbone atom in m that will be added to a Chain
   for (i = 0; i < m.natoms(); i++) {
      if (-1 != index_map[i]) {
         m_end = i;
         break;
      }
   }

   for (i = 0; i < m.natoms(); i++) {
      if (-1 == index_map[i])
         continue;
      ti = m.type(i);

      // Interactions of particle i with other particles in m
      for (j = i+1; j < m.natoms(); j++) {
         if (-1 == index_map[j])
            continue;
         if (m.bond_separation(i, j) <= _bonded_cutoff)
            continue;  // i, j "bonded", interact via torsion energy

         tj = m.type(j);
         r = m.position[j] - m.position[i];  // unwrapped but in same monomer
         rsq = r.sqlength();
         if (rsq < energy.pair_range_sq())
            E += energy.pair_energy(ti, tj, rsq);
      }

      // Now we start to consider interactions with particles already in Cells,
      // but we don't want to include pair interactions with particles that 
      // are within _bonded_cutoff bonds; these interactions were included in
      // the torsion selection.  Any particle connected to the current tail 
      // particle by cell_bonds or fewer should NOT interact with particle i
      cell_bonds = _bonded_cutoff - m.bond_separation(m_end, i) - 1;
      if (cell_bonds < 0)
         cell_bonds =  0;

      // Map particle i to a Cell in sys
      pos = m.position[i];
      sys.wrap(pos);
      icell = sys.cell_index(pos);

      // Interactions of particle i with particles in Cell icell
      for (j = 0; j < sys.cells[icell].nparticles(); j++) {
         j_ndx = sys.cells[icell].index(j);
         if ((cell_bonds > 0) && 
             (ic == sys.cells[icell].chain(j)) && 
             (itail > -1) &&
             (sys.bond_separation(c, itail, j_ndx, cell_bonds))) {
            continue;  // particles i and j interact only via torsion energy
         }
         tj = sys.cells[icell].ptype(j);
         r = sys.cells[icell].pos[j] - pos;  // wrapped
         rsq = r.sqlength();
         if (rsq < energy.pair_range_sq())
            E += energy.pair_energy(ti, tj, rsq);
      }

      // Interactions of particle i with particles in Cells surrounding icell
      sys.find_neighbor_cells(icell, nbr_cells);
      for (k = 0; k < 26; k++) {
         if (-1 == nbr_cells[k])
            continue;
         for (j = 0; j < sys.cells[nbr_cells[k]].nparticles(); j++) {
            j_ndx = sys.cells[nbr_cells[k]].index(j);
            if ((cell_bonds > 0) && 
                (ic == sys.cells[nbr_cells[k]].chain(j)) && 
                (itail > -1) &&
                (sys.bond_separation(c, itail, j_ndx, cell_bonds))) {
               continue;  // particles i and j interact only via torsion energy
            }
            tj = sys.cells[nbr_cells[k]].ptype(j);
            rsq = sys.min_sqdist(pos, sys.cells[nbr_cells[k]].pos[j]); //wrapped
            if (rsq < energy.pair_range_sq())
               E += energy.pair_energy(ti, tj, rsq);
         }
      }

      // Interactions with excluded regions
      for (iexc = 0; iexc < exclusions.size(); iexc++) {
         if (!exclusions[iexc])
            continue;
         if (exclusions[iexc]->inside(pos))
            E += BIG_VAL;
         else {
            rsq = exclusions[iexc]->min_sqdist(pos);
            if (rsq < energy.pair_range_sq())
               E += energy.hard_energy(ti, ti, rsq);
         }
      }
   }
   return E;
}

// ============================================================================
// Return an energy for the indicated torsion angle from linear interpolation 
// between _torsion_energy[] values
// ============================================================================
Real
CBMC::_interpolate_energy(Real angle) const
{
   Real dang, dE;
   const int nang = _torsion_angle.size()-1;
   int i;

   for (i = 0; i < nang; i++) {
      if ((angle >= _torsion_angle[i]) && (angle < _torsion_angle[i+1]))
         break;
   }
   if (i == nang) {
      dang = _torsion_angle[0] + 360.0 - _torsion_angle[i];
      dE = _torsion_energy[0] - _torsion_energy[i];
   }
   else {
      dang = _torsion_angle[i+1] - _torsion_angle[i];
      dE = _torsion_energy[i+1] - _torsion_energy[i];
   }
   return _torsion_energy[i] + (angle - _torsion_angle[i])*dE/dang;
}

// ============================================================================
// Set torsion angles on the backbone of the Monomer m to minimize interactions
// with surrounding particles
// ============================================================================
void
CBMC::_rotate_torsions(Monomer& m, 
                       const vector<int> index_map, 
                       Chain& c, 
                       int ic, 
                       int itail, 
                       bool first, 
                       System& sys, 
                       const vector<const Exclusion *>& exclusions,
                       const Energy& energy,
                       RNG& rng) const
{
   Real Ef, Ei, dE, Enb;
   int i, j, k;
   const int nbb = m.nbackbone();
   // The initial torsion to rotate along the backbone of m: if this is the
   // first monomer, then we start with 3, so we can use 4 actual particle
   // positions (0, 1, 2, 3); otherwise, start with 2 (previous tail, 0, 1, 2)
   int j_init = (first) ? 3 : 2;  
   int j_final = nbb;
   vector<Real> prev_angle(nbb);

   if (_only_between) {
      // torsion nbb-2 controls the rotation about the bond between monomers
      j_init = nbb-2;
      j_final = nbb-1;
      // if first monomer, no rotation in this case
      if (first)
         j_init = j_final;
   }

   Ei = 0.0;
   for (j = j_init; j < j_final; j++) {
      k = rng.select_by_weight(_torsion_prob);
      m.set_torsion(j, _torsion_angle[k]);
      Ei += _torsion_energy[k];
   }
   m.update_positions();
   Enb = _nonbond_interactions(m, index_map, c, ic, itail, sys, exclusions, 
                               energy);
   if (_log_energies)
      log_message("   Initial energy: torsion %g, pair %g\n", Ei, Enb);
   if (Enb > BIG_VAL)
      log_message("ERROR: initial configuration in exclusion!\n");
   Ei += Enb;

   for (i = 0; i < _nconfigs; i++) {
      Ef = 0.0;
      for (j = j_init; j < j_final; j++) {
         prev_angle[j] = m.torsion_angle(j);
         k = rng.select_by_weight(_torsion_prob);
         m.set_torsion(j, _torsion_angle[k]);
         Ef += _torsion_energy[k];
      }
      m.update_positions();
      Enb = _nonbond_interactions(m, index_map, c, ic, itail, sys, exclusions, 
                                  energy);
      Ef += Enb;
      dE = Ef - Ei;
      if ((Enb > BIG_VAL) || ((dE > 0.0) && (exp(-dE*_inv_kT) < rng.yield()))) {
         // Reject this configuration
         for (j = j_init; j < j_final; j++)
            m.set_torsion(j, prev_angle[j]);
      }
      else {  // Keep this configuration
         if (_log_energies) {
            log_message("   Configuration %d energy: torsion %g, pair %g\n", 
                        i+1, Ef-Enb, Enb);
         }
         Ei = Ef;
      }
   }

   for (i = 0; i < _ndelta_steps; i++) {
      Ef = 0.0;
      for (j = j_init; j < j_final; j++) {
         prev_angle[j] = m.torsion_angle(j);
         if (rng.yield() < 0.5)
            m.set_torsion(j, prev_angle[j] + _delta_degs);
         else
            m.set_torsion(j, prev_angle[j] - _delta_degs);
         Ef += _interpolate_energy(m.torsion_angle(j));
      }
      m.update_positions();
      Enb = _nonbond_interactions(m, index_map, c, ic, itail, sys, exclusions, 
                                  energy);
      Ef += Enb;
      dE = Ef - Ei;
      if ((Enb > BIG_VAL) || ((dE > 0.0) && (exp(-dE*_inv_kT) < rng.yield()))) {
         // Reject this configuration
         for (j = j_init; j < j_final; j++)
            m.set_torsion(j, prev_angle[j]);
      }
      else  { // Keep this configuration
         if (_log_energies) {
            log_message("   Step %d energy: torsion %g, pair %g\n", 
                        i+1, Ef-Enb, Enb);
         }
         Ei = Ef;
      }
   }

   // Final positions
   m.update_positions();
}

// ============================================================================
// Prepare Chain data members
// ============================================================================
void
CBMC::_resize_chain_data()
{
   _pattern_index.reserve(_nchains);
   _chain_tail_anchor.reserve(_nchains);
   _chain_tail_pos.reserve(_nchains);
   _chain_tail_vec.reserve(_nchains);
   for (int i = 0; i < _nchains; i++) {
      _pattern_index.push_back(0);
      _chain_tail_anchor.push_back(Triple());
      _chain_tail_pos.push_back(Triple());
      _chain_tail_vec.push_back(Triple());
   }
}

// ============================================================================
// Add a new bond; add to BOTH particles
// ============================================================================
void
CBMC::_add_bond(System& sys,
                Chain& ch,
                Ivec& bond) const
{
   int cell, slot;

   //printf("adding bond %d - %d\n", bond[0], bond[1]); fflush(stdout);
   cell = ch.cell(bond[0]);
   slot = ch.slot(bond[0]);
   sys.cells[cell].add_bond(slot, bond[1], bond[2]);
   cell = ch.cell(bond[1]);
   slot = ch.slot(bond[1]);
   sys.cells[cell].add_bond(slot, bond[0], bond[2]);
}

// ============================================================================
// Build chains 
//
// This is the main Chain building method for this Builder.
// ============================================================================
void 
CBMC::build(System& sys,
            vector<Chain>& chains, 
            vector<Monomer>& monomers,
            vector<int>& monomer_selection,
            const vector<AMOCH::Chaintype>& chain_types, 
            const vector<const Exclusion *>& exclusions,
            const Energy& energy,
            RNG& rng)
{
   int i, ic, im, chtype, i0;
   int nadded = 0;
   int nbuild;
   int nbb, Nch, nmon;
   int cell, slot;
   int offset;
   int nmonpart;
   int new_tail_index;
   unsigned iexc;
   bool first, last;
   Real len, theta, cosphi, g;
   Real maxsqlen;
   Triple pos, v;
   const Triple size = sys.size();
   const Triple min = sys.min();
   Ivec bond(3);

   nmon = 0;
   len = 0.0;
   for (im = 0; im < (int)monomers.size(); im++) {
      if (monomers[im].natoms() > nmon)
         nmon = monomers[im].natoms();
      if (monomers[im].length() > len)
         len = monomers[im].length();
   }
   maxsqlen = len*len;  // max square length of any monomer type

   vector<int> index_map(nmon, -1);

   if (0 == _pattern_index.size()) {  // if not restarted
      _nchains = chains.size();
      _resize_chain_data();
      //
      // Initial chain positions, orientations
      //
      for (ic = 0; ic < _nchains; ic++) {
         // Randomly position first particle in Chain ic; should be at least
         // a monomer length away from any excluded region
         do {
            _chain_tail_pos[ic].x = min.x + rng.yield() * size.x;
            _chain_tail_pos[ic].y = min.y + rng.yield() * size.y;
            _chain_tail_pos[ic].z = min.z + rng.yield() * size.z;
            for (iexc = 0; iexc < exclusions.size(); iexc++) {
               if (exclusions[iexc]->min_sqdist(_chain_tail_pos[ic]) < maxsqlen)
                  break;
            }
         }  while (iexc < exclusions.size());

         // By convention, next particle is oriented along the X axis wrt the
         // first particle; here we rotate away from X to a random angle
         theta = 2.0 * rng.yield() * M_PI;
         cosphi = 2.0 * rng.yield() - 1.0;
         g = _bb_bondlen * sqrt(1.0 - cosphi*cosphi);
         v.x = g * cos(theta);
         v.y = g * sin(theta);
         v.z = _bb_bondlen * cosphi;
         _chain_tail_vec[ic] = v;
         _chain_tail_vec[ic].normalize();

         // Again by convention, the anchor position for the torsion calculation
         // associated with the third particle is oriented along Y; here we 
         // apply the same rotation to Y
         _chain_tail_anchor[ic].x = v.y;
         g = 1.0/(v.x + 1.0);
         _chain_tail_anchor[ic].y = -v.x - (v.z*v.z)*g;
         _chain_tail_anchor[ic].z = v.y*v.z*g;
      }
   }

   if (_nbuild > -1)
      nbuild = _nbuild;
   else {
      nbuild = 0;
      for (ic = 0; ic < _nchains; ic++) 
         nbuild += chains[ic].max_monomers();
   }

   //
   // Main build loop over all monomers added
   //
   while (nadded < nbuild) {

      // Select chain
      ic = -1;
      while (-1 == ic) {
         ic = (int)(rng.yield() * _nchains);
         if (chains[ic].nmonomers() == chains[ic].max_monomers()) {
            log_message("Searching for chain to build; ");
            log_message("chain %06d is complete\n", ic+1);
            ic = -1;
         }
      }

      first = (chains[ic].nmonomers() == 0);
      last = (chains[ic].nmonomers() == chains[ic].max_monomers()-1);
      chtype = chains[ic].type();

      // Select monomer
      im = -1;
      if (last) 
         im = chain_types[chtype].term_monomer_index();
      if (im == -1) {  // no terminating Monomer
         if (chain_types[chtype].select_by_weight()) 
            im = chain_types[chtype].monomer_index_by_weight(rng);
         else {  // pattern
            if (_pattern_index[ic] == chain_types[chtype].nmonomers())
               _pattern_index[ic] = 0;
            im = chain_types[chtype].monomer_index(_pattern_index[ic]++);
         }
      }
      monomer_selection[im]++;
      nbb = monomers[im].nbackbone();

      //
      // Add monomers[im] to chains[ic]
      //

      // Number of particles already in Chain ic
      Nch = chains[ic].nparticles();

      // Number of particles from Monomer im to be added to Chain ic
      nmon = (_bb_only) ? nbb : monomers[im].natoms();

      // Map indices of atoms in Monomer to Chain particles
      for (i = 0; i < (int)index_map.size(); i++)
         index_map[i] = -1;

      nmonpart = 0;
      offset = Nch - 1;
      if (first) {
         index_map[0] = 0;
         nmonpart++;
         offset++;
      }
      for (i = 1; i < nbb-1; i++) {
         index_map[i] = offset + i;
         nmonpart++;
      }
      // i == nbb-1
      new_tail_index = index_map[i-1];
      if (last) {
         index_map[i] = offset + i;
         new_tail_index = index_map[i];
         nmonpart++;
      }
      else 
         offset--;

      for (i = nbb; i < nmon; i++) {
         index_map[i] = offset + i;
         nmonpart++;
      }

      // Set (unwrapped) positions in monomers[im]
      monomers[im].set_anchor(_chain_tail_anchor[ic]);
      monomers[im].position[0] = _chain_tail_pos[ic];
      len = (first) ? monomers[im].bond_length(nbb-1) : _bb_bondlen;
      _chain_tail_vec[ic] *= len;
      monomers[im].position[1] = monomers[im].position[0] + _chain_tail_vec[ic];

      if (_rotate) {
         _rotate_torsions(monomers[im], index_map, chains[ic], ic, 
                          chains[ic].tail_index, first, sys, exclusions, 
                          energy, rng);
      }
      else
         monomers[im].update_positions();

      // Record selected torsion values along backbone; if this is the first
      // monomer, don't record the first torsion (2) -- with only 3 atoms, 
      // any rotation will involve an artifical point.  Better to accept the
      // torsion defined in the monomer input.
      i0 = (first) ? 3 : 2;
      for (i = i0; i < nbb; i++)
         chains[ic].record_torsion(monomers[im].torsion_angle(i));

      // Record unwrapped tail positions
      _chain_tail_anchor[ic] = monomers[im].position[nbb-3];
      _chain_tail_pos[ic] = monomers[im].position[nbb-2];
      _chain_tail_vec[ic] = 
         monomers[im].position[nbb-1] - monomers[im].position[nbb-2];
      _chain_tail_vec[ic].normalize();

      // Add new particles to system (with empty bonds)
      for (i = 0; i < nmon; i++) {
         if (-1 == index_map[i])
            continue;
         // Monomer positions are unwrapped
         pos = monomers[im].position[i];
         sys.wrap(pos);
         cell = sys.cell_index(pos);
         slot = sys.cells[cell].add_particle(ic, index_map[i], 
                                             monomers[im].type(i),
                                             pos,
                                             monomers[im].position[i]);
         chains[ic].add_particle_cell(index_map[i], cell, slot);
      }
      // Add bond between monomers
      offset = nmon - 3;
      if (!first) {
         bond[0] = chains[ic].tail_index;
         bond[1] = index_map[1];
         bond[2] = _bb_bondtype;
         _add_bond(sys, chains[ic], bond);
      }
      // Add new bonds to system
      for (i = 0; i < monomers[im].nbonds(); i++) {
         bond = monomers[im].bond(i);  // [index1, index2, type]
         if ((bond[0] >= nmon) || (bond[1] >= nmon))
            continue;  // bond i not used in system
         if ((index_map[bond[0]] == -1) || (index_map[bond[1]] == -1))
            continue;  // bond i not used in system
         /*
         if ((bond[0] == 0) && !first) {
            bond[1] = index_map[bond[1]];
            bond[0] = bond[1] - offset;
         }
         else if ((bond[1] == 0) && !first) {
            bond[0] = index_map[bond[0]];
            bond[1] = bond[0] - offset;
         }
         else if ((bond[0] == nbb-1) && !last) {
            bond[1] = index_map[bond[1]];
            bond[0] = bond[1] + offset;
         }
         else if ((bond[1] == nbb-1) && !last) {
            bond[0] = index_map[bond[1]];
            bond[1] = bond[0] + offset;
         }
         else {
            bond[0] = index_map[bond[0]];
            bond[1] = index_map[bond[1]];
         }
         */
         bond[0] = index_map[bond[0]];
         bond[1] = index_map[bond[1]];
         _add_bond(sys, chains[ic], bond);
      }
      chains[ic].add_monomer(nmonpart);
      chains[ic].tail_index = new_tail_index;
      nadded++;
      if (_log_monomers) {
         log_message("Add monomer %s to chain %d; %d / %d in chain, ", 
                     monomers[im].name().c_str(), ic+1, chains[ic].nmonomers(),
                     chains[ic].max_monomers());
         log_message("%d / %d overall\n", nadded, nbuild);
      }
      if (_log_status)
         AMOCH::show_status(AMOCH::get_logfile(), nadded, nbuild);
   }  // end main build loop
}

// ============================================================================
// Write internal state to the output FILE f
// ============================================================================
void 
CBMC::write(FILE *f) const
{
   fwrite((void *)&_nchains, sizeof(int), 1, f);
   fwrite((void *)&_pattern_index[0], sizeof(int), _nchains, f);
   fwrite((void *)&_chain_tail_anchor[0], sizeof(Triple), _nchains, f);
   fwrite((void *)&_chain_tail_pos[0], sizeof(Triple), _nchains, f);
   fwrite((void *)&_chain_tail_vec[0], sizeof(Triple), _nchains, f);
}

// ============================================================================
// Read internal state from the input FILE f
// ============================================================================
void 
CBMC::read(FILE *f)
{
   fread((void *)&_nchains, sizeof(int), 1, f);
   _resize_chain_data();
   fread((void *)&_pattern_index[0], sizeof(int), _nchains, f);
   fread((void *)&_chain_tail_anchor[0], sizeof(Triple), _nchains, f);
   fread((void *)&_chain_tail_pos[0], sizeof(Triple), _nchains, f);
   fread((void *)&_chain_tail_vec[0], sizeof(Triple), _nchains, f);
}
