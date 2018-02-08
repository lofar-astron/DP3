# $Id$

include(LofarPackageVersion)
include(PythonInstall)

lofar_add_library(dppp
  Package__Version.cc
  DPRun.cc DPStep.cc DPInput.cc DPBuffer.cc DPInfo.cc
  DPLogger.cc ProgressMeter.cc FlagCounter.cc
  UVWCalculator/UVWCalculator.cc  BaselineSelection.cc ApplyCal.cc
  MSReader.cc MultiMSReader.cc MSWriter.cc MSUpdater.cc Counter.cc
  Averager.cc MedFlagger.cc PreFlagger.cc UVWFlagger.cc
  StationAdder.cc ScaleData.cc Filter.cc PhaseShift.cc
  Demixer.cc 
  Position.cc Stokes.cc SourceDBUtil.cc
  Apply.cc EstimateMixed.cc EstimateNew.cc Simulate.cc Simulator.cc
  SubtractMixed.cc SubtractNew.cc
  ModelComponent.cc PointSource.cc GaussianSource.cc Patch.cc
  ModelComponentVisitor.cc GainCal.cc StefCal.cc
  DemixerNew.cc DemixInfo.cc DemixWorker.cc
  Predict.cc OneApplyCal.cc
  ApplyBeam.cc
  PhaseFitter.cc H5Parm.cc SolTab.cc DummyStep.cc H5ParmPredict.cc
  Upsample.cc
)

lofar_add_bin_program(NDPPP NDPPP.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)

install(CODE "EXECUTE_PROCESS(
  COMMAND \"${CMAKE_COMMAND}\" -E create_symlink NDPPP DPPP WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin
  )"
)

lofar_add_bin_scripts(taqlflagger)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)
