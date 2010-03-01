# $Id$

include(LofarPackageVersion)

lofar_add_library(dppp
  Package__Version.cc
  DPRun.cc DPStep.cc DPInput.cc DPBuffer.cc AverageInfo.cc
  ProgressMeter.cc FlagCounter.cc
  UVWCalculator.cc
  MSReader.cc MSWriter.cc MSUpdater.cc
  Averager.cc MedFlagger.cc PreFlagger.cc
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
  Pipeline.cc)

lofar_add_bin_program(NDPPP DPPP.cc)
lofar_add_bin_program(IDPPP IDPPP.cc PipelineProcessControl.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)

install(PROGRAMS
  mwflagger
  mwflagger-part
  taqlflagger
  DESTINATION bin)
