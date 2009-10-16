# $Id$

include(LofarPackageVersion)

lofar_add_library(dppp
  Package__Version.cc
  FlaggerStatistics.cc
  MsInfo.cc
  MsFile.cc
  ComplexMedianFlagger.cc
  FrequencyFlagger.cc
  ComplexMedianFlagger2.cc
  MADFlagger.cc
  BandpassCorrector.cc
  DataSquasher.cc
  DataBuffer.cc
  TimeBuffer.cc
  RunDetails.cc
  ProgressMeter.cc
  Pipeline.cc)

lofar_add_bin_program(CS1_IDPPP IDPPP.cc PipelineProcessControl.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)

install(PROGRAMS
  mwflagger
  mwflagger-part
  DESTINATION bin)
