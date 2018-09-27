cmake_minimum_required(VERSION 3.0)

project(DP3)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

if(NOT DISABLE_OMP)
  find_package( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(NOT DISABLE_OMP)

find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIRS})

find_package(Casacore COMPONENTS casa ms tables REQUIRED)
include_directories(${CASACORE_INCLUDE_DIR})
# Add old casacore include directory because LOFAR beam library did
# not use the casacore/ prefix in the includes until 3.2 (Sep '18)
include_directories(${CASACORE_INCLUDE_DIR}/casacore)

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/Station.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
  include_directories(${LOFAR_STATION_RESPONSE_DIR})
  add_definitions(-DHAVE_LOFAR_BEAM)
  message("LOFAR beam library found.")
  
  set(LOFAR_DEPENDENT_FILES
    DPPP/ApplyBeam.cc
    DPPP/DemixerNew.cc 
    DPPP/DemixInfo.cc
    DPPP/DemixWorker.cc 
  )
  
else()
  if(WITH_LOFAR_BEAM)
    message(FATAL_ERROR "LOFAR beam library not found")
  else()
    message("LOFAR beam library not found -- disabling LOFAR beam prediction.")
    set(LOFAR_STATION_RESPONSE_LIB "")
    set(LOFAR_DEPENDENT_FILES "")
  endif(WITH_LOFAR_BEAM)
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Object libraries (which are just file groups, not actual libraries) are made in order
# to avoid compiling everything multiple times for multiple executables.
add_library(Common OBJECT 
  Common/BaselineSelect.cc
  Common/ClusterDesc.cc
  Common/DataConvert.cc
  Common/NodeDesc.cc
  Common/ParameterHandler.cc
  Common/ParameterRecord.cc
  Common/ParameterSet.cc
  Common/ParameterSetImpl.cc
  Common/ParameterValue.cc
  Common/PrettyUnits.cc
  Common/StringUtil.cc
  Common/Timer.cc
  Common/TypeNames.cc
  Common/VdsDesc.cc
  Common/VdsMaker.cc
  Common/VdsPartDesc.cc
)
set(COMMON_OBJECT $<TARGET_OBJECTS:Common>)

add_library(ParmDB OBJECT
  ParmDB/Axis.cc
  ParmDB/AxisMapping.cc
  ParmDB/Box.cc
  ParmDB/Grid.cc
  ParmDB/Parm.cc
  ParmDB/ParmCache.cc
  ParmDB/ParmDB.cc
  ParmDB/ParmDBBlob.cc
  ParmDB/ParmDBCasa.cc
  ParmDB/ParmDBLocker.cc
  ParmDB/ParmDBMeta.cc
  ParmDB/ParmFacade.cc
  ParmDB/ParmFacadeLocal.cc
  ParmDB/ParmFacadeRep.cc
  ParmDB/ParmSet.cc
  ParmDB/ParmValue.cc
  ParmDB/PatchInfo.cc
  ParmDB/SourceData.cc
  ParmDB/SourceDB.cc
  ParmDB/SourceDBBlob.cc
  ParmDB/SourceDBCasa.cc
  ParmDB/SourceInfo.cc
)
set(PARMDB_OBJECT $<TARGET_OBJECTS:ParmDB>)

add_library(Blob OBJECT
  Blob/BlobAipsIO.cc
  Blob/BlobArray.cc
  Blob/BlobHeader.cc
  Blob/BlobIBufStream.cc
  Blob/BlobIStream.cc
  Blob/BlobOBufStream.cc
  Blob/BlobOStream.cc
)
set(BLOB_OBJECT $<TARGET_OBJECTS:Blob>)

add_library(DPPP_OBJ OBJECT
  DPPP/DPRun.cc DPPP/DPStep.cc DPPP/DPInput.cc DPPP/DPBuffer.cc
  DPPP/DPInfo.cc DPPP/DPLogger.cc DPPP/ProgressMeter.cc DPPP/FlagCounter.cc
  DPPP/UVWCalculator.cc  DPPP/BaselineSelection.cc DPPP/ApplyCal.cc
  DPPP/MSReader.cc DPPP/MultiMSReader.cc DPPP/MSWriter.cc DPPP/MSUpdater.cc
  DPPP/Counter.cc DPPP/Averager.cc DPPP/MedFlagger.cc DPPP/PreFlagger.cc
  DPPP/UVWFlagger.cc DPPP/StationAdder.cc DPPP/ScaleData.cc DPPP/Filter.cc 
  DPPP/PhaseShift.cc DPPP/Demixer.cc DPPP/Position.cc DPPP/Stokes.cc 
  DPPP/SourceDBUtil.cc DPPP/Apply.cc DPPP/EstimateMixed.cc DPPP/EstimateNew.cc 
  DPPP/Simulate.cc DPPP/Simulator.cc DPPP/SubtractMixed.cc DPPP/SubtractNew.cc
  DPPP/ModelComponent.cc DPPP/PointSource.cc DPPP/GaussianSource.cc DPPP/Patch.cc
  DPPP/ModelComponentVisitor.cc DPPP/GainCal.cc DPPP/StefCal.cc
  DPPP/Predict.cc DPPP/OneApplyCal.cc
  DPPP/PhaseFitter.cc DPPP/H5Parm.cc DPPP/SolTab.cc 
  DPPP/DummyStep.cc DPPP/H5ParmPredict.cc DPPP/GridInterpolate.cc DPPP/Upsample.cc
  DPPP/Split.cc
  ${LOFAR_DEPENDENT_FILES}
)
set(DPPP_OBJECT $<TARGET_OBJECTS:DPPP_OBJ>)

add_subdirectory(DPPP)
add_subdirectory(TestDynDPPP)
add_subdirectory(PythonDPPP)
add_subdirectory(DPPP_AOFlag)
# SPW Combine has dependencies to LOFAR
#add_subdirectory(SPWCombine)
add_subdirectory(DPPP_Interpolate)

find_package(Armadillo)
if(${ARMADILLO_FOUND})
  add_subdirectory(DPPP_DDECal)
else()
  message(WARNING "Armadillo was not found, NOT building DPPP_DDECal")
endif()

# Required for the ObjectFactory
include_directories(Common)

add_executable(makesourcedb ParmDB/makesourcedb.cc ${PARMDB_OBJECT} ${BLOB_OBJECT} ${COMMON_OBJECT})
target_link_libraries(makesourcedb ${CASACORE_LIBRARIES} ${Boost_SYSTEM_LIBRARY})

install (TARGETS makesourcedb DESTINATION bin) 
