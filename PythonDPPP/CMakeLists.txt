# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# We set Boost-python as required in toplevel build file.
if(TRUE)

  include(PythonInstall)

  add_library(_pythondppp MODULE PythonDPPP.cc)
  set_target_properties(_pythondppp PROPERTIES
    SOVERSION ${DPPP_GLOBAL_SOVERSION}
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_BUILD_DIR}/lofar/pythondppp)

# This is a quick-and-dirty fix to install the Python binding module in the
# right place. It will now be installed twice, because lofar_add_library()
# will install it in $prefix/$libdir
  install(TARGETS _pythondppp
    DESTINATION ${PYTHON_INSTALL_DIR}/lofar/pythondppp)

# We only need this library because dependent packages need to link against
# this package's version info.
#add_library(dppp_pythondppp
#  PythonStep.cc
#  DPStepBase.cc
#)

  # Install Python modules
  python_install(__init__.py DESTINATION lofar/pythondppp)

  find_package(Casacore COMPONENTS python3)

#  add_subdirectory(test)
endif()
