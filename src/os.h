// ============================================================================
// os.h -- Prototypes for OS-specific functions defined in <os>.cc
// ----------------------------------------------------------------------------
// Copyright (c) 2019 Benjamin P. Haley 
// ----------------------------------------------------------------------------
// See the LICENSE file for information on usage and redistribution of this
// file and for a DISCLAIMER OF ALL WARRANTIES. 
// ============================================================================

#ifndef AMOCH_OS_H
#define AMOCH_OS_H

#include <string>

namespace AMOCH {

// Return a joined path e.g. "dir/file"
std::string join_path(const std::string& dir,
                      const std::string& file);

// Return true if path exists as a file, else return false
bool
check_file(const std::string& path);

// Return true if dir exists as a directory, else return false
bool
check_dir(const std::string& dir);

// Create a new directory; return true on success, false on error after logging
// a message
bool create_dir(const std::string& dir);

}   // AMOCH

#endif  // AMOCH_OS_H
