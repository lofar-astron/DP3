cmake_minimum_required(VERSION 2.8.8)

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

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/Station.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
	include_directories(${LOFAR_STATION_RESPONSE_DIR})
  set(LBEAM_LIBS ${LOFAR_STATION_RESPONSE_LIB})
  set(LBEAM_FILES lofar/lbeamevaluator.cpp lofar/lmspredicter.cpp aterms/lofarbeamterm.cpp)
  message("LOFAR beam library found.")
else()
  message(FATAL_ERROR "LOFAR beam library not found")
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

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
