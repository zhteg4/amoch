// ============================================================================
// monomer.cc -- AMOCH::Monomer methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include <cstdlib>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include "monomer.h"
#include "logfile.h"
#include "defs.h"

using std::string;
using std::vector;
using OpenBabel::OBAtom;
using OpenBabel::OBConversion;
using OpenBabel::OBInternalCoord;
using OpenBabel::OBBondIterator;
using AMOCH::Monomer;
using AMOCH::Real;
using AMOCH::Triple;
using AMOCH::Ivec;
using AMOCH::PartType;
using AMOCH::BondType;
using AMOCH::ZEntry;
using AMOCH::log_message;
using AMOCH::Param;

// ============================================================================
// ZEntry constructor
// ============================================================================
ZEntry::ZEntry(int typ,
               int bat,
               Real blen,
               int aat,
               Real bang,
               int tat,
               Real tang) :
   type(typ),
   bond_atom(bat),
   bond_length(blen),
   angle_atom(aat),
   bond_angle(bang),
   torsion_atom(tat),
   torsion_angle(tang)
{}

// ============================================================================
// Monomer onstructor
// ============================================================================
Monomer::Monomer(const string& mname) :
   _name(mname),
   _mol(),
   _natoms(0),
   _nbb(0),
   _nbonds(0),
   _anchor(),
   _mass(0.0),
   _bb_mass(0.0),
   _length(0.0),
   _bonds(0, Ivec()),
   _zmatrix(0, ZEntry()),
   position(0, Triple())
{}

// ============================================================================
// Read molecule from file.  Return true on success, false on error with 
// message written to log file.
// ============================================================================
bool
Monomer::read(const char *path)
{
   OBConversion conv;

   if (!conv.SetInFormat(OBConversion::FormatFromExt(path))) {
      log_message("Unknown input format for molecule file %s\n", path);
      return false;
   }

   std::ifstream ifs(path);

   if (!ifs.is_open()) {
      log_message("Unable to open molecule file %s\n", path);
      return false;
   }

   bool status = conv.Read(&_mol, &ifs);

   ifs.close();
   if (!status) {
      log_message("Error reading monomer file %s\n", path);
      return false;
   }
   
   _natoms = _mol.NumAtoms();
   log_message("Read monomer %s from file %s\n", _name.c_str(), path);
   log_message("   %d atoms in molecule\n", _natoms);

   return status;
}

// ============================================================================
// Add non-backbone atoms to id and reorder molecule
// ============================================================================
void
Monomer::_reorder(vector<int> id)
{
   int i, j, n;

   // Now we have the backbone; add remaining atoms to id
   n = _nbb;
   for (i = 0; i < _natoms; i++) {
      for (j = 0; j < n; j++) {
         if (id[j] == i)
            break;
      }
      if (j == n)  // atom i not in backbone
         id[n++] = i;
   }
   // Reorder atoms: atom ids (start at 0) in id[] need to be converted to
   // indices (start at 1)
   for (i = 0; i < _natoms; i++)
      id[i]++;
   _mol.RenumberAtoms(id); // reassigns idx but not id; use GetIdx() later...
}

// ============================================================================
// Identify backbone atoms from head to tail;  reorder molecule so backbone 
// atoms are first, followed by all side group atoms.
// ============================================================================
void
Monomer::identify_backbone(int head,
                           int tail)
{
   vector<OBBondIterator> bit(_natoms);
   vector<int> id(_natoms);
   OBAtom *bbat, *nbrat;
   int i;
   bool have_bb = false;

   // Head and tail are 0-based indices; internally we use atom id, which 
   // also starts at 0, from OBAtom::GetId(), to identify the backbone 
   _nbb = 1;
   id[0] = head;
   bbat = _mol.GetAtomById(head);
   nbrat = bbat->BeginNbrAtom(bit[0]);
   while (!have_bb) {
      while ((nbrat)) {
         for (i = 0; i < _nbb; i++) {
            if (id[i] == (int)nbrat->GetId())
               break;
         }
         if (i == _nbb) {  // nbrat is NOT in backbone; add it
            id[_nbb] = (int)nbrat->GetId();
            bbat = nbrat;
            nbrat = bbat->BeginNbrAtom(bit[_nbb++]);
            break;
         }
         else  // nbrat is already in backbone
            nbrat = bbat->NextNbrAtom(bit[_nbb-1]);
      }
      if ((int)bbat->GetId() == tail)
         have_bb = true;
      if (!nbrat) {
         --_nbb;
         // All neighbors of bbat exhausted; reject bbat and move on to next
         // neighbor of previous backbone atom
         bbat = _mol.GetAtomById(id[_nbb-1]);
         nbrat = bbat->NextNbrAtom(bit[_nbb-1]);
      }
   }
   log_message("   %d backbone atoms\n", _nbb);
   _reorder(id);
}

// ============================================================================
// Specify backbone atoms from head to tail;  reorder molecule so backbone 
// atoms are first, followed by all side group atoms.
// ============================================================================
void
Monomer::specify_backbone(const vector<int>& backbone)
{
   vector<int> id(_natoms);

   for (_nbb = 0; _nbb < (int)backbone.size(); _nbb++)
      id[_nbb] = backbone[_nbb]-1;
   log_message("   %d backbone atoms\n", _nbb);
   _reorder(id);
}

// ============================================================================
// Set iternal coordinates, along with positions, mass, and length; 
// add new particle types to part_types.  Return true on success, or false
// on error after logging a message.
// ============================================================================
bool
Monomer::setup(const Param& pmon,
               vector<PartType>& part_types)
{
   OBAtom *a;
   Triple d;
   Real len;
   int i, n;

   // Positions, mass
   position.reserve(_natoms);
   for (i = 0; i < _natoms; i++) {
      a = _mol.GetAtom(i+1);
      position.push_back(Triple(a->GetX(), a->GetY(), a->GetZ()));
      _mass += a->GetAtomicMass();
      if (i < _nbb)
         _bb_mass += a->GetAtomicMass();
      for (int j = 0; j < i; j++) {
         d = position[i] - position[j];
         len = d.length();
         if (len > _length)
            _length = len;
      }
   }
   log_message("   Molecule mass: %g amu\n", _mass);

   // Internal coordinates
   int el, nb, hb, res, typ;
   int bat, aat, tat;

   vector<OBInternalCoord *> vic(_natoms+1, NULL);
   for (i = 0; i < _natoms; i++)
      vic[i+1] = new OBInternalCoord;
   OpenBabel::CartesianToInternal(vic, _mol);

   _zmatrix.reserve(_natoms);
   for (i = 0; i < _natoms; i++) {
      n = i+1;
      a = _mol.GetAtom(n);
      el = a->GetAtomicNum();
      nb = a->GetValence();
      hb = a->IsHbondDonorH();
      res = a->IsInRing();
      typ = AMOCH::match_part_type(part_types, el, nb, hb, res);
      if (-1 == typ) {
         typ = part_types.size();
         part_types.push_back(PartType(el, OpenBabel::etab.GetSymbol(el), nb, 
                              hb, res));
      }
      bat = (vic[n]->_a) ? vic[n]->_a->GetIdx()-1 : -1;
      aat = (vic[n]->_b) ? vic[n]->_b->GetIdx()-1 : -1;
      tat = (vic[n]->_c) ? vic[n]->_c->GetIdx()-1 : -1;
      _zmatrix.push_back(ZEntry(typ, bat, vic[n]->_dst,
                                     aat, vic[n]->_ang,
                                     tat, vic[n]->_tor));
   }
   if (_natoms > 2)
      _zmatrix[2].torsion_angle = 90.0;  // OB convention, if not specified
   for (i = 0; i < _natoms; i++)
      delete vic[i+1];

   FILE *f = AMOCH::get_logfile();

   if ((f)) {
      log_message("Internal coordinates:\n");
      write_zmatrix(f, part_types);
   }

   // set torsions {...}
   const Param *p;

   if ((p = pmon.find("set"))) {
      if (p->value() == "torsions") {
         const Param *pang = p->find("angle");
         const Param *pind = p->find("index");
         Real ang;

         if (!pang) {
            log_message("\nMissing \"angle\" values for \"set torsions\"\n");
            return false;
         }
         if (!pind) {
            log_message("\nMissing \"index\" values for \"set torsions\"\n");
            return false;
         }
         if (pind->nvalues() != pang->nvalues()) {
            log_message("\n\"angle\" and \"index\" have different sizes ");
            log_message("for \"set torsions\"\n");
            return false;
         }
         for (i = 0; i < pind->nvalues(); i++) {
            n = atoi(pind->value(i).c_str()) - 1;
            ang = atof(pang->value(i).c_str());
            _zmatrix[n].torsion_angle = ang;
            log_message("Set torsion %d to %g degrees\n", n+1, ang);
         }
         if ((f)) {
            log_message("Updated internal coordinates:\n");
            write_zmatrix(f, part_types);
         }
      }
      else {
         log_message("\nUnknown \"set\" option, %s, in monomer %s\n", 
                        p->value().c_str(), _name.c_str());
         return false;
      }
   }
   return true;
}

// ============================================================================
// Identify bonds; if internal coordinates are set, report extra bonds not 
// in _zmatrix.  Add new bond types to bond_types
// ============================================================================
void 
Monomer::find_bonds(vector<BondType>& bond_types)
{
   OpenBabel::OBBond *ob;
   OBBondIterator bit;
   int i, btype, ptype1, ptype2, order;

   if (0 == _mol.NumBonds())
      _mol.ConnectTheDots();  // let OB identify bonds
   _nbonds = _mol.NumBonds();

   _bonds.reserve(_nbonds);
   for (i = 0; i < _nbonds; i++)
      _bonds.push_back(Ivec(3));

   ob = _mol.BeginBond(bit);
   i = 0;
   while ((ob)) {
      _bonds[i][0] = ob->GetBeginAtomIdx()-1;
      ptype1 = _zmatrix[_bonds[i][0]].type;
      _bonds[i][1] = ob->GetEndAtomIdx()-1;
      ptype2 = _zmatrix[_bonds[i][1]].type;
      order = (int)ob->GetBondOrder();
      btype = AMOCH::match_bond_type(bond_types, ptype1, ptype2, order);
      if (-1 == btype) {
         btype = bond_types.size();
         bond_types.push_back(AMOCH::BondType(ptype1, ptype2, order));
      }
      _bonds[i][2] = btype;
      ob = _mol.NextBond(bit);
      i++;
   }
   if ((_nbonds > _natoms-1) && (_zmatrix.size() > 0)) {
      log_message("Extra bonds, not listed in z-matrix:\n");
      for (i = 0; i < _nbonds; i++) {
         if ((_zmatrix[_bonds[i][0]].bond_atom != _bonds[i][1]) && 
             (_zmatrix[_bonds[i][1]].bond_atom != _bonds[i][0]))
            log_message("   %d - %d\n", _bonds[i][0]+1, _bonds[i][1]+1);
      }
   }
}

// ============================================================================
// Adjust internal coordinates to ensure torsion rotations preserve structure
// ============================================================================
void
Monomer::adjust_internal_coords(const vector<PartType>& part_types)
{
   if (0 == _zmatrix.size())
      return;

   int i, j, k, ba;
   int prev_angle_atom;
   Triple bai, bak;
   Triple n1, n2;

   for (i = _nbb; i < _natoms; i++) {
      ba = _zmatrix[i].bond_atom;
      // ba is bonded to i
      for (j = 0; j < _nbonds; j++) {
         k = -1;
         if (_bonds[j][0] == ba)
            k = _bonds[j][1];
         else if (_bonds[j][1] == ba)
            k = _bonds[j][0];
         // ba is bonded to k ...
         if ((k != i) && (k != _zmatrix[i].torsion_atom) && 
             (k != _zmatrix[i].angle_atom) && (k > 1) && (k < _nbb)) {
            // ...and k can rotate: make i angle depend on k
            prev_angle_atom = _zmatrix[i].angle_atom;
            _zmatrix[i].angle_atom = k;
            bai = position[i] - position[ba];
            bai.normalize();
            bak = position[k] - position[ba];
            bak.normalize();
            _zmatrix[i].bond_angle = RAD2DEG(acos(dot(bai, bak)));
            // XXX hack?  search for new torsion_atom?
            _zmatrix[i].torsion_atom = prev_angle_atom;
            n1 = cross(bai, bak);
            n1.normalize();
            bai = position[prev_angle_atom] - position[ba];
            n2 = cross(bai, bak);
            n2.normalize();
            _zmatrix[i].torsion_angle = RAD2DEG(acos(dot(n1, n2)));
            // XXX hack?
            if ((i > _nbb) && 
     (fabs(_zmatrix[i].torsion_angle-_zmatrix[i-1].torsion_angle) < 1.0)) {
               _zmatrix[i].torsion_angle *= -1.0;
            }
            break;
         }
      }
   }

   FILE *f = AMOCH::get_logfile();

   if ((f)) {
      log_message("Adjusted internal coordinates:\n");
      write_zmatrix(f, part_types);
   }
}

// ============================================================================
// Calculate the position of an atom given the internal coordinates
//
// This function implements the position calculation of OpenBabel's 
// InternalToCartesian(). 
// ============================================================================
static Triple
calc_pos(const Triple& a,
         const Triple& b,
         const Triple& c,
         Real blen,
         Real bang,
         Real tang)
{
   Triple v1, v2, v3, n, nn;

   v1 = a - b;
   v2 = a - c;
   n = cross(v1, v2);
   nn = cross(v1, n);
   n.normalize();
   nn.normalize();
   n *= -sin(tang);
   nn *= cos(tang);
   v3 = n + nn;
   v3.normalize();
   v3 *= blen*sin(bang);
   v1.normalize();
   v1 *= blen*cos(bang);
   return a + v3 - v1;
}

// ============================================================================
// Update atomic positions from internal coordinates; if tail_length > 0.0,
// use it as the length of the tail atom bond.  
//
// Requires position[0], position[1], and _anchor
// ============================================================================
void 
Monomer::update_positions(Real tail_length)
{
   Real blen, bang, tang;
   Triple a, b, c;
   Real tblen = _zmatrix[_nbb-1].bond_length;

   if (tail_length > 0.0)
      _zmatrix[_nbb-1].bond_length = tail_length;

   // Positions of atoms 0 and 1 are already set; for atom 2 we use the
   // anchor position for the torsion calculation.
   a = position[_zmatrix[2].bond_atom] - position[0];
   blen = _zmatrix[2].bond_length;
   b = position[_zmatrix[2].angle_atom] - position[0];
   bang = DEG2RAD(_zmatrix[2].bond_angle);
   c = _anchor - position[0];
   //tang = DEG2RAD(90.0);
   tang = DEG2RAD(_zmatrix[2].torsion_angle);
   position[2] = calc_pos(a, b, c, blen, bang, tang) + position[0];
   // For all other atoms we have a torsion value
   for (int i = 3; i < _natoms; i++) {
      a = position[_zmatrix[i].bond_atom] - position[0];
      blen = _zmatrix[i].bond_length;
      b = position[_zmatrix[i].angle_atom] - position[0];
      bang = DEG2RAD(_zmatrix[i].bond_angle);
      c = position[_zmatrix[i].torsion_atom] - position[0];
      tang = DEG2RAD(_zmatrix[i].torsion_angle);
      position[i] = calc_pos(a, b, c, blen, bang, tang) + position[0];
   }
   _zmatrix[_nbb-1].bond_length = tblen;
}

// ============================================================================
// Return the number of bonds between atoms m and n (m < n); returns _nbonds
// if no path is found (error)
// ============================================================================
int 
Monomer::bond_separation(int m,
                         int n) const
{
   vector<int> bonded(_natoms);
   int mbonds = 0;
   int nbonds = 0;
   int i = _zmatrix[n].bond_atom;
   int j;

   if (m == n)
      return 0;
   while (i != -1) {
      bonded[nbonds++] = i;
      if (i == m)
         return nbonds;
      i = _zmatrix[i].bond_atom;
   }
   i = _zmatrix[m].bond_atom;
   while (i != -1) {
      mbonds++;
      for (j = 0; j < nbonds; j++) {
         if (i == bonded[j])
            return mbonds+j+1;
      }
      i = _zmatrix[i].bond_atom;
   }
   return _nbonds;  // FIXME ...
}

// ============================================================================
// Return the maximum number of bonds on any atom in the Monomer m
// ============================================================================
int 
Monomer::max_bonds_per_atom() const
{
   int bmax = 0;
   int b;

   for (int i = 0; i < _natoms; i++) {
      b = 0;
      for (int j = 0; j < _nbonds; j++) {
         if ((_bonds[j][0] == i) || (_bonds[j][1] == i))
            b++;
      }
      if (b > bmax)
         bmax = b;
   }
   return bmax;
}

// ============================================================================
// Write the zmatrix to the output FILE f
// ============================================================================
void 
Monomer::write_zmatrix(FILE *f,
                       const vector<PartType>& part_types) const
{
   int i;

   for (i = 0; i < _natoms; i++) {
      fprintf(f, "%3.3s", part_types[_zmatrix[i].type].symbol().c_str());
      if (_zmatrix[i].bond_atom > -1) {
         fprintf(f, "%4d%6.2f", _zmatrix[i].bond_atom+1, 
                 _zmatrix[i].bond_length);
      }
      if (_zmatrix[i].angle_atom > -1) {
         fprintf(f, "%4d%9.2f", _zmatrix[i].angle_atom+1, 
                 _zmatrix[i].bond_angle);
      }
      if (_zmatrix[i].torsion_atom > -1) {
         fprintf(f, "%4d%9.2f", _zmatrix[i].torsion_atom+1, 
                 _zmatrix[i].torsion_angle);
      }
      fprintf(f, "\n");
   }
   fflush(f);
}

// ============================================================================
// Write atomic postitions to an output file
// ============================================================================
void
Monomer::write(const char *path)
{
   OBConversion conv;

   if (!conv.SetOutFormat(OBConversion::FormatFromExt(path))) {
      log_message("Unknown input format for monomer dump file %s\n", path);
      return;
   }

   std::ofstream ofs(path);

   if (!ofs.is_open()) {
      log_message("Unable to open monomer dump file %s\n", path);
      return;
   }

   OBAtom *a;

   for (int i = 0; i < _natoms; i++) {
      a = _mol.GetAtom(i+1);
      a->SetVector(position[i].x, position[i].y, position[i].z);
   }
   (void)conv.Write(&_mol, &ofs);
   ofs.close();
}

