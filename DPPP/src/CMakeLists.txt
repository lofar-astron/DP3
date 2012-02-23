# $Id$

include(LofarPackageVersion)
include(PythonInstall)

lofar_add_library(dppp
  Package__Version.cc
  DPRun.cc DPStep.cc DPInput.cc DPBuffer.cc DPInfo.cc
  ParSet.cc DPLogger.cc ProgressMeter.cc FlagCounter.cc
  UVWCalculator.cc
  MSReader.cc MultiMSReader.cc MSWriter.cc MSUpdater.cc Counter.cc
  Averager.cc MedFlagger.cc AORFlagger.cc PreFlagger.cc UVWFlagger.cc
  PhaseShift.cc Demixer.cc BBSExpr.cc EstimateNDPPP.cc
  FlaggerStatistics.cc
# The following is old IDPPP stuff.
#  MsInfo.cc
#  MsFile.cc
#  ComplexMedianFlagger.cc
#  FrequencyFlagger.cc
#  ComplexMedianFlagger2.cc
#  MADFlagger.cc
#  BandpassCorrector.cc
#  DataSquasher.cc
#  DataBuffer.cc
#  TimeBuffer.cc
#  RunDetails.cc
#  Pipeline.cc
)

lofar_add_bin_program(NDPPP NDPPP.cc)
#lofar_add_bin_program(IDPPP IDPPP.cc PipelineProcessControl.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)

install(PROGRAMS
  mwflagger
  mwflagger-part
  taqlflagger
  DESTINATION bin)

# Python modules.
python_install(
  __init__.py 
  DESTINATION lofar/dppp)
