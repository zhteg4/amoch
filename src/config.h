// ============================================================================
// src/config.h.  Generated from config.h.in by configure.
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_CONFIG_H
#define AMOCH_CONFIG_H

#include <cfloat>

#define REAL_DOUBLE

namespace AMOCH {

// Floating point type
#ifdef REAL_DOUBLE

typedef double Real;
#define REAL_SMALL         DBL_EPSILON*1e4
#define REAL_MAX           DBL_MAX
#define SCN_REAL           "%lf"
#define BIG_VAL            1e70

#else

typedef float Real;
#define REAL_SMALL         FLT_EPSILON*1e2
#define REAL_MAX           FLT_MAX
#define SCN_REAL           "%f"
#define BIG_VAL            1e35

#endif

}   // AMOCH

// Mersenne exponent controls the period of RNG sequences
#define SFMT_MEXP          19937

// Cell size is this much larger than pair interaction range
#define RANGE_EPS          0.0

// Standard file names
#define ATOMS_DATA_FILE    "atoms.dat"
#define ATOM_TYPE_FILE     "atom_type.dat.in"
#define BONDS_DATA_FILE    "bonds.dat"
#define BOND_TYPE_FILE     "bond_type.dat.in"
#define CHAIN_DATA_FILE    "chain%06d.dat"
#define RNG_DATA_FILE      "rng.dat"
#define LAMMPS_DUMP_FILE   "polymer.dump"
#define BUILDER_DATA_FILE  "builder.dat"

// Scale factor used to resize vectors
#define RESIZE_SCALE       1.3

#endif  // AMOCH_CONFIG_H
