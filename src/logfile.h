// ============================================================================
// logfile.h -- Prototypes for message logging functions defined in logfile.cc
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#ifndef AMOCH_LOGFILE_H
#define AMOCH_LOGFILE_H

#include <cstdio>

namespace AMOCH {

// Set log file to f
void set_logfile(FILE *f);

// Open log file
void open_logfile(const char *path);

// Return log file pointer
FILE *get_logfile(void);

// Write message to log file
void log_message(const char *msg,
                 ...);

// Close log file, if open
void close_logfile(void);

}   // AMOCH

#endif  // AMOCH_LOGFILE_H
