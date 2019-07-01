// ============================================================================
// os.cc -- OS-specific functions 
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#include <cstring>
#include <cerrno>
#include "os.h"
#include "logfile.h"

using std::string;
using AMOCH::log_message;

#ifdef __unix

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

// ============================================================================
// Return a joined path e.g. "dir/file"
// ============================================================================
string 
AMOCH::join_path(const string& dir,
                 const string& file)
{
   return dir + "/" + file;
}

// ============================================================================
// Return true if path exists as a file, else return false
// ============================================================================
bool
AMOCH::check_file(const std::string& path)
{
   FILE *f = fopen(path.c_str(), "r");

   if ((f)) {
      fclose(f);
      return true;
   }
   return false;
}

// ============================================================================
// Return true if dir exists, else return false
// ============================================================================
bool
AMOCH::check_dir(const string& dir)
{
   bool ret = false;
   const char *cdir = dir.c_str();
   DIR *d = opendir(cdir);

   if ((d)) {
      ret = true;
      closedir(d);
   }
   else {
      if (errno != ENOENT)
         log_message("Unable to check %s: %s\n", cdir, strerror(errno));
   }
   return ret;
}

// ============================================================================
// Create a new directory; return true on success, false on error after logging
// a message
// ============================================================================
bool
AMOCH::create_dir(const string& dir)
{
   const char *cdir = dir.c_str();

   if (0 != mkdir(cdir, 0755)) {
      log_message("Unable to create directory %s: %s\n", cdir, strerror(errno));
      return false;
   }
   return true;
}
#endif  // __unix
