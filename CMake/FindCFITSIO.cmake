# - Try to find CFITSIO.
# Variables used by this module:
#  CFITSIO_ROOT_DIR     - CFITSIO root directory
# Variables defined by this module:
#  CFITSIO_FOUND        - system has CFITSIO
#  CFITSIO_INCLUDE_DIR  - the CFITSIO include directory (cached)
#  CFITSIO_INCLUDE_DIRS - the CFITSIO include directories
#                         (identical to CFITSIO_INCLUDE_DIR)
#  CFITSIO_LIBRARY      - the CFITSIO library (cached)
#  CFITSIO_LIBRARIES    - the CFITSIO libraries
#                         (identical to CFITSIO_LIBRARY)

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

if(NOT CFITSIO_FOUND)

  if(NOT CFITSIO_ROOT_DIR AND DEFINED ENV{CFITSIO_ROOT_DIR})
    set(CFITSIO_ROOT_DIR "$ENV{CFITSIO_ROOT_DIR}")
  endif()

  find_path(CFITSIO_INCLUDE_DIR fitsio.h
    HINTS ${CFITSIO_ROOT_DIR} PATH_SUFFIXES include include/cfitsio include/libcfitsio0)
  find_library(CFITSIO_LIBRARY cfitsio
    HINTS ${CFITSIO_ROOT_DIR} PATH_SUFFIXES lib)
  find_library(M_LIBRARY m)
  mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY M_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(CFITSIO DEFAULT_MSG
    CFITSIO_LIBRARY M_LIBRARY CFITSIO_INCLUDE_DIR)

  set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY} ${M_LIBRARY})

endif(NOT CFITSIO_FOUND)
