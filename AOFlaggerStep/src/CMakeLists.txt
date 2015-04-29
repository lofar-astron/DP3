# $Id: CMakeLists.txt 30439 2014-11-19 15:04:34Z dijkema $

include(LofarPackageVersion)

lofar_add_library(aoflaggerstep
  Package__Version.cc
  AOFlaggerStep.cc
)

lofar_add_bin_program(versionaoflaggerstep versionaoflaggerstep.cc)

#target_link_libraries (aoflaggerstep dppp aoflagger)
