# $Id: CMakeLists.txt 30439 2014-11-19 15:04:34Z dijkema $

include(LofarPackageVersion)
include(PythonInstall)

lofar_add_library(_pythondppp MODULE PythonDPPP.cc)
set_target_properties(_pythondppp PROPERTIES 
  PREFIX ""
  LIBRARY_OUTPUT_DIRECTORY ${PYTHON_BUILD_DIR}/lofar/pythondppp)

# This is a quick-and-dirty fix to install the Python binding module in the
# right place. It will now be installed twice, because lofar_add_library()
# will install it in $prefix/$libdir
install(TARGETS _pythondppp
  DESTINATION ${PYTHON_INSTALL_DIR}/lofar/pythondppp)

# We only need this library because dependent packages need to link against
# this package's version info.
lofar_add_library(dppp_pythondppp
  Package__Version.cc
  PythonStep.cc
  DPStepBase.cc
)

lofar_add_bin_program(versionpythondppp versionpythondppp.cc)

# Install Python modules
python_install(__init__.py DESTINATION lofar/pythondppp)
