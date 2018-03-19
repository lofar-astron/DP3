# $Id: CMakeLists.txt 30439 2014-11-19 15:04:34Z dijkema $

include(LofarPackageVersion)

lofar_add_library(dppp_sagecal
  Package__Version.cc
  SageCal.cc Register.cc
)

lofar_add_bin_program(versiondppp_sagecal versiondppp_sagecal.cc)
