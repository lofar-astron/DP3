# FindHDF5 uses NATIVE_COMMAND in separate_arguments, which requires
# CMake 3.9.
cmake_minimum_required(VERSION 3.9)

project(DP3)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3")

if(PORTABLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native") 
endif(PORTABLE)

enable_testing()

find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIRS})

find_package(Casacore COMPONENTS casa ms tables fits REQUIRED)
include_directories(${CASACORE_INCLUDE_DIR})
# Add old casacore include directory because LOFAR beam library did
# not use the casacore/ prefix in the includes until 3.2 (Sep '18)
include_directories(${CASACORE_INCLUDE_DIR}/casacore)

find_package(Threads REQUIRED)

find_package(IDGAPI NO_MODULE QUIET)

if(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
#  set(IDG_FILES ....)
  include_directories(${IDGAPI_INCLUDE_DIRS})
  set(HAVE_IDG TRUE)
  add_definitions(-DHAVE_IDG)
  message(STATUS "Image domain gridder API libraries found.")
else(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  set(IDGAPI_LIBRARIES)
  set(IDG_FILES)
  message(WARNING "Image domain gridder API libraries NOT found. IDG prediction will not be available.")
endif(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)

find_package(CFITSIO REQUIRED)

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time filesystem python program_options system unit_test_framework REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/Station.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
  include_directories(${LOFAR_STATION_RESPONSE_DIR})
  add_definitions(-DHAVE_LOFAR_BEAM)
  message(STATUS "LOFAR beam library found.")
  
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

add_library(IDGPredict OBJECT
  IDGPredict/FitsIOChecker.cc
  IDGPredict/FitsReader.cc
  IDGPredict/FitsWriter.cc)
  
set(IDGPREDICT_OBJECT $<TARGET_OBJECTS:IDGPredict>)

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
  DPPP/ModelComponentVisitor.cc DPPP/GainCal.cc DPPP/GainCalAlgorithm.cc
  DPPP/Predict.cc DPPP/OneApplyCal.cc
  DPPP/PhaseFitter.cc DPPP/H5Parm.cc DPPP/SolTab.cc 
  DPPP/DummyStep.cc DPPP/H5ParmPredict.cc DPPP/GridInterpolate.cc DPPP/Upsample.cc
  DPPP/Split.cc
  DPPP/Interpolate.cc
  ${LOFAR_DEPENDENT_FILES}
)
set(DPPP_OBJECT $<TARGET_OBJECTS:DPPP_OBJ>)

#
# Dependencies of AOFlagger
#
find_package(PythonInterp REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(PNG REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIR})

find_library(AOFLAGGER_LIB aoflagger REQUIRED)
find_path(AOFLAGGER_INCLUDE_DIR NAMES aoflagger.h)
include_directories(${AOFLAGGER_INCLUDE_DIR})

find_package(PkgConfig)
pkg_check_modules(GTKMM gtkmm-3.0>=3.0.0)
pkg_check_modules(SIGCXX sigc++-2.0)
if(GTKMM_FOUND)
  set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} ${GTKMM_LIBRARIES} ${GLIBMM_LIBRARIES})
endif(GTKMM_FOUND)

add_library(AOFlaggerStep_OBJ OBJECT
  AOFlaggerStep/AOFlaggerStep.cc
  AOFlaggerStep/Register.cc
)
set(AOFLAGGERSTEP_OBJECT $<TARGET_OBJECTS:AOFlaggerStep_OBJ>)

set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES}
  ${AOFLAGGER_LIB}
  ${CASACORE_LIBRARIES}
  ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}  ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_PYTHON_LIBRARY}
  ${EXTRA_LIBRARIES}
  ${PYTHON_LIBRARIES}
  ${IDGAPI_LIBRARIES}
  Threads::Threads)
  
#
# DDECal dependencies
#
find_package(Armadillo)
if(${ARMADILLO_FOUND})
  add_definitions(-DHAVE_ARMADILLO)
  set(DDE_ARMADILLO_FILES
    DDECal/KLFitter.cc 
    DDECal/PiercePoint.cc 
    DDECal/ScreenConstraint.cc)
else()
  message(WARNING "Armadillo was not found, not including screenfitter inside DDECal")
  set(DDE_ARMADILLO_FILES)
endif()

include_directories(${ARMADILLO_INCLUDE_DIRS})
set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} ${ARMADILLO_LIBRARIES})
add_library(DDECal_OBJ OBJECT
  DDECal/DDECal.cc DDECal/Register.cc DDECal/Stopwatch.cc
  DDECal/MultiDirSolver.cc
  DDECal/Constraint.cc
  DDECal/SmoothnessConstraint.cc
  DDECal/TECConstraint.cc DDECal/RotationConstraint.cc
  DDECal/RotationAndDiagonalConstraint.cc
  ${DDE_ARMADILLO_FILES})
set(DDECAL_OBJECT $<TARGET_OBJECTS:DDECal_OBJ>)

set(DPPP_GLOBAL_SOVERSION 2)


add_subdirectory(DPPP)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  add_subdirectory(PythonDPPP)
endif()
# SPW Combine has dependencies to LOFAR
#add_subdirectory(SPWCombine)

# Required for the ObjectFactory
include_directories(Common)

add_executable(makesourcedb ParmDB/makesourcedb.cc ${PARMDB_OBJECT} ${BLOB_OBJECT} ${COMMON_OBJECT})
target_link_libraries(makesourcedb ${CASACORE_LIBRARIES} ${Boost_SYSTEM_LIBRARY})

add_executable(showsourcedb ParmDB/showsourcedb.cc ${PARMDB_OBJECT} ${BLOB_OBJECT} ${COMMON_OBJECT})
target_link_libraries(showsourcedb ${CASACORE_LIBRARIES} ${Boost_SYSTEM_LIBRARY})

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  SET(OS_SPECIFIC_TESTS) # No specific tests for Apple
else()
  # These run only on Linux
  SET(OS_SPECIFIC_TESTS 
#    TestDynDPPP/test/tDynStep.cc # This test still fails
  )
endif()

add_library(dppp_testdyndppp SHARED
  TestDynDPPP/TestDynStep.cc ${DPPP_OBJECT} ${PARMDB_OBJECT} ${BLOB_OBJECT} ${DDECAL_OBJECT} ${COMMON_OBJECT} ${AOFLAGGERSTEP_OBJECT}
)
set_target_properties(dppp_testdyndppp PROPERTIES VERSION ${DPPP_GLOBAL_SOVERSION})

# I'm avoiding cmake's test commands, because they build executable tests by default
# which makes the compilation process more expensive.
set(TEST_FILENAMES
  DPPP/test/tMirror.cc
  DPPP/test/tPhaseShift.cc
  DDECal/test/tRotationConstraint.cc
)

set(RUNTEST_SOURCES
  DPPP/test/runtests.cc
  ${TEST_FILENAMES}
  ${OS_SPECIFIC_TESTS}
  ${DPPP_OBJECT} ${COMMON_OBJECT} ${BLOB_OBJECT} ${DDECAL_OBJECT} ${IDGPREDICT_OBJECT} ${PARMDB_OBJECT} ${AOFLAGGERSTEP_OBJECT}
)

add_executable(runtests EXCLUDE_FROM_ALL ${RUNTEST_SOURCES})

target_link_libraries(runtests ${EXTRA_LIBRARIES} ${LOFAR_STATION_RESPONSE_LIB})

add_test(buildruntests "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target runtests)
# Add each test as a cmake test (only for tests added with BOOST_AUTO_TEST_SUITE and BOOST_AUTO_TEST_CASE)
foreach(TEST_FILENAME ${TEST_FILENAMES})
  file(READ ${TEST_FILENAME} TEST_FILE_CONTENTS)

  # Extract test suite name
  string(REGEX MATCHALL "BOOST_AUTO_TEST_CASE\\( *([A-Za-z_0-9]+) *\\)" FOUND_TESTS ${TEST_FILE_CONTENTS})
  string(REGEX MATCH "BOOST_AUTO_TEST_SUITE\\(*([A-Za-z_0-9]+)*\\)" TESTSUITE_NAME_FULL ${TEST_FILE_CONTENTS})

  # Extract test name
  string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+) *\\).*" "\\1" TESTSUITE_NAME ${TESTSUITE_NAME_FULL})
  foreach(FOUND_TEST ${FOUND_TESTS})
    string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+) *\\).*" "\\1" TEST_NAME ${FOUND_TEST})
    add_test(NAME "${TESTSUITE_NAME}.${TEST_NAME}" COMMAND runtests --run_test=${TESTSUITE_NAME}/${TEST_NAME} --catch_system_error=yes)
    set_tests_properties(${TESTSUITE_NAME}.${TEST_NAME} PROPERTIES DEPENDS buildruntests)
  endforeach()
endforeach()

# I'll leave this here as an example how separate test executables could be added
#add_executable(tMirror EXCLUDE_FROM_ALL
#  DPPP/test/tMirror.cc
#)
#target_link_libraries(tMirror ${CASACORE_LIBRARIES})

#add_executable(tPhaseShift EXCLUDE_FROM_ALL
#  DPPP/test/tPhaseShift.cc
#  ${DPPP_OBJECT} ${COMMON_OBJECT} ${BLOB_OBJECT} ${DDECAL_OBJECT} ${PARMDB_OBJECT} ${AOFLAGGERSTEP_OBJECT})
#target_link_libraries(tPhaseShift ${EXTRA_LIBRARIES} ${LOFAR_STATION_RESPONSE_LIB})

install (TARGETS makesourcedb showsourcedb DESTINATION bin) 
