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

find_library (CDDLIB_LIBRARY NAMES cdd)

set (CDDLIB_LIBRARIES ${CDDLIB_LIBRARY})
set (CDDLIB_INCLUDE_DIRS ${CDDLIB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cddlib DEFAULT_MSG CDDLIB_LIBRARY CDDLIB_INCLUDE_DIRS)

mark_as_advanced (
  CDDLIB_LIBRARY
  CDDLIB_LIBRARIES
  CDDLIB_INCLUDE_DIR
  CDDLIB_INCLUDE_DIRS)

