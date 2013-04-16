# $Id$

include(LofarPackageVersion)
include(PythonInstall)

lofar_add_library(dppp
  Package__Version.cc
  DPRun.cc DPStep.cc DPInput.cc DPBuffer.cc DPInfo.cc
  DPLogger.cc ProgressMeter.cc FlagCounter.cc
  UVWCalculator.cc  BaselineSelection.cc
  MSReader.cc MultiMSReader.cc MSWriter.cc MSUpdater.cc Counter.cc
  Averager.cc MedFlagger.cc AORFlagger.cc PreFlagger.cc UVWFlagger.cc
  StationAdder.cc ScaleData.cc Filter.cc
  PhaseShift.cc Demixer.cc
  Position.cc Stokes.cc SourceDBUtil.cc
  Apply.cc EstimateMixed.cc Simulate.cc Simulator.cc SubtractMixed.cc
  ModelComponent.cc PointSource.cc GaussianSource.cc Patch.cc
  ModelComponentVisitor.cc
  FlaggerStatistics.cc
)

lofar_add_bin_program(NDPPP NDPPP.cc)
#lofar_add_bin_program(IDPPP IDPPP.cc PipelineProcessControl.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)
lofar_add_bin_program(measures measures.cc)

install(PROGRAMS
  mwflagger
  mwflagger-part
  taqlflagger
  DESTINATION bin)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)
