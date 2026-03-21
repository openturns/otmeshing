# - Find cddlib
# A general dimension code for computing convex hulls and related structures
# http://www.qhull.org/
#
# The module defines the following variables:
#  CDDLIB_INCLUDE_DIRS, where to find libqhull_r/qhull_ra.h, etc.
#  CDDLIB_LIBRARIES, the libraries needed to use cddlib.
#  CDDLIB_FOUND, If false, do not try to use cddlib.
# also defined, but not for general use are
#  CDDLIB_LIBRARY, where to find the cddlib library.
#
#=============================================================================
# Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path (CDDLIB_INCLUDE_DIR cdd.h PATH_SUFFIXES cddlib)

# version
set (_VERSION_FILE ${CDDLIB_INCLUDE_DIR}/cddtypes.h)
if (EXISTS ${_VERSION_FILE})
  #define ddf_DDVERSION   "Version 0.94m"
  file (STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "#define dd_DDVERSION")
  if (_VERSION_LINE)
    string (REGEX REPLACE ".*\"Version (.*)\"" "\\1" cddlib_VERSION ${_VERSION_LINE})
  endif ()
endif ()

# check version
set (_CDDLIB_VERSION_MATCH TRUE)
if (cddlib_FIND_VERSION AND cddlib_VERSION)
  if (cddlib_FIND_VERSION_EXACT)
    if (${cddlib_FIND_VERSION} VERSION_EQUAL ${cddlib_VERSION})
    else()
      set (_CDDLIB_VERSION_MATCH FALSE)
    endif ()
  else ()
    if (${cddlib_FIND_VERSION} VERSION_GREATER ${cddlib_VERSION})
      set (_CDDLIB_VERSION_MATCH FALSE)
    endif ()
  endif ()
endif ()

find_library (CDDLIB_LIBRARY NAMES cdd)

set (CDDLIB_LIBRARIES ${CDDLIB_LIBRARY})
set (CDDLIB_INCLUDE_DIRS ${CDDLIB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cddlib
  REQUIRED_VARS CDDLIB_LIBRARY CDDLIB_INCLUDE_DIRS
  VERSION_VAR cddlib_VERSION)

mark_as_advanced (
  CDDLIB_LIBRARY
  CDDLIB_LIBRARIES
  CDDLIB_INCLUDE_DIR
  CDDLIB_INCLUDE_DIRS
  cddlib_VERSION)

