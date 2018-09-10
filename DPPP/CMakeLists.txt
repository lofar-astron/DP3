include(PythonInstall)

# Required for the ObjectFactory
include_directories(../Common)

add_library(dppp
  DPRun.cc DPStep.cc DPInput.cc DPBuffer.cc DPInfo.cc
  DPLogger.cc ProgressMeter.cc FlagCounter.cc
  UVWCalculator.cc  BaselineSelection.cc ApplyCal.cc
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
  PhaseFitter.cc H5Parm.cc SolTab.cc DummyStep.cc H5ParmPredict.cc GridInterpolate.cc
  Upsample.cc Split.cc
  
  ../Blob/BlobAipsIO.cc
  ../Blob/BlobArray.cc
  ../Blob/BlobHeader.cc
  ../Blob/BlobIBufStream.cc
  ../Blob/BlobIStream.cc
  ../Blob/BlobOBufStream.cc
  ../Blob/BlobOStream.cc
  
  ../Common/BaselineSelect.cc
  ../Common/ClusterDesc.cc
  ../Common/DataConvert.cc
  ../Common/NodeDesc.cc
  ../Common/ParameterHandler.cc
  ../Common/ParameterRecord.cc
  ../Common/ParameterSet.cc
  ../Common/ParameterSetImpl.cc
  ../Common/ParameterValue.cc
  ../Common/PrettyUnits.cc
  ../Common/StringUtil.cc
  ../Common/Timer.cc
  ../Common/TypeNames.cc
  ../Common/VdsDesc.cc
  ../Common/VdsMaker.cc
  ../Common/VdsPartDesc.cc
  
  ../ParmDB/Axis.cc
  ../ParmDB/AxisMapping.cc
  ../ParmDB/Box.cc
  ../ParmDB/Grid.cc
  ../ParmDB/Parm.cc
  ../ParmDB/ParmCache.cc
  ../ParmDB/ParmDB.cc
  ../ParmDB/ParmDBBlob.cc
  ../ParmDB/ParmDBCasa.cc
  ../ParmDB/ParmDBLocker.cc
  ../ParmDB/ParmDBMeta.cc
  ../ParmDB/ParmFacade.cc
  ../ParmDB/ParmFacadeLocal.cc
  ../ParmDB/ParmFacadeRep.cc
  ../ParmDB/ParmSet.cc
  ../ParmDB/ParmValue.cc
  ../ParmDB/PatchInfo.cc
  ../ParmDB/SourceData.cc
  ../ParmDB/SourceDB.cc
  ../ParmDB/SourceDBBlob.cc
  ../ParmDB/SourceDBCasa.cc
  ../ParmDB/SourceInfo.cc
)
target_link_libraries(dppp)

add_executable(DPPP NDPPP.cc)
target_link_libraries(DPPP dppp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${CASACORE_LIBRARIES} ${LOFAR_STATION_RESPONSE_LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# lofar_add_bin_scripts(taqlflagger)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)

add_subdirectory(share)
#add_subdirectory(test)
