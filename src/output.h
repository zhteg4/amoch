// ============================================================================
// output.h -- Prototypes for functions defined in output.cc
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_OUTPUT_H
#define AMOCH_OUTPUT_H

#include "system.h"
#include "types.h"

namespace AMOCH {

// Write a status bar with completion percentage to the output FILE f
void show_status(FILE *f,
                 int step,
                 int nsteps);
      
// Write a PDB output file of all Chains in a System
void write_pdb(FILE *f, 
               const char *comment,
               const System& sys,
               const std::vector<Chain>& chains,
               const std::vector<PartType>& part_types,
               bool wrapped,
               bool connect = true);

// Write a XYZ output file of all Chains in a System
void write_xyz(FILE *f, 
               const char *comment,
               const System& sys,
               const std::vector<Chain>& chains,
               const std::vector<PartType>& part_types,
               bool wrapped);

// Write atoms.dat file for data4Lammps
void 
write_atoms_dat(FILE *f,
                const System& sys, 
                const std::vector<Chain>& chains);
                 
// Write atom_type.dat input file for data4Lammps; the particle types are not
// specific to any force field, but all the info is written for subsequent
// forcefield typing
void 
write_atom_type_dat_in(FILE *f, 
                       const std::vector<PartType>& part_types);

// Write bonds.dat input file for data4Lammps
void
write_bonds_dat(FILE *f, 
                const System& sys, 
                const std::vector<Chain>& chains);

// Write bond_type.dat input file for data4Lammps; the atom types are written
// as indices for types in atom_type.dat, not specific to any force field
void 
write_bond_type_dat_in(FILE *f,
                       const std::vector<AMOCH::BondType>& bond_types);

}   // AMOCH

#endif  // AMOCH_OUTPUT_H
