# $Id: CMakeLists.txt 30439 2014-11-19 15:04:34Z dijkema $

include(LofarPackageVersion)

lofar_add_library(dppp_aoflag
  Package__Version.cc
  AOFlaggerStep.cc
#  SlidingFlagger.cc
  Register.cc
)

lofar_add_bin_program(versiondppp_aoflag versiondppp_aoflag.cc)
