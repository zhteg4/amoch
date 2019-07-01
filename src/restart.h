// ============================================================================
// restart.h -- Prototypes for functions defined in restart.cc
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_RESTART_H
#define AMOCH_RESTART_H

#include "system.h"
#include "chain.h"

namespace AMOCH {

// Return true if particles in LAMMPS dump file specified by path can be read
// into sys, else log an error message and return false.  This function expects
// a specific column order and UNwrapped coordinates; the dump file can be 
// written using these arguments with a custom dump: "id mol type xu yu zu"
// For example,
//    dump 1 all custom 1 file.dump id mol type xu yu zu
bool read_lammps_dump(FILE *f,
                      System& sys,
                      std::vector<Chain>& chains);

// Return true if the bond info in the indicated file can be read into sys,
// else return false after logging a message
bool read_bonds(FILE *f,
                System& sys,
                std::vector<Chain>& chains);

}   // AMOCH

#endif  // AMOCH_RESTART_H
