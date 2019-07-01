// ============================================================================
// logfile.cc -- Message logging functions
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES.
// ============================================================================

#include <cstdarg>
#include "logfile.h"

// Log file
static FILE *flog = NULL;

// ============================================================================
// Set log file to f
// ============================================================================
void
AMOCH::set_logfile(FILE *f)
{
   flog = f;
}

// ============================================================================
// Open log file
// ============================================================================
void 
AMOCH::open_logfile(const char *path)
{
   FILE *f = fopen(path, "w");

   if ((f))
      flog = f;
}

// ============================================================================
// Return log file pointer
// ============================================================================
FILE *
AMOCH::get_logfile(void)
{
   return flog;
}

// ============================================================================
// Write message to log file
// ============================================================================
void 
AMOCH::log_message(const char *msg,
                   ...)
{
   if ((flog) && (msg)) {
      va_list ap;

      va_start(ap, msg);
      vfprintf(flog, msg, ap);
      va_end(ap);
      fflush(flog);
   }
}

// ============================================================================
// Close log file, if open
// ============================================================================
void 
AMOCH::close_logfile(void)
{
   if ((flog) && (flog != stdout) && (flog != stderr)) {
      fclose(flog);
      flog = NULL;
   }
}

