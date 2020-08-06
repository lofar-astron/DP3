# FindHDF5 uses NATIVE_COMMAND in separate_arguments, which requires
# CMake 3.9.
cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0074 NEW)

project(DP3)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

add_compile_options(-Wall -O3)

if(NOT PORTABLE)
  add_compile_options(-march=native)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include(CTest)
enable_testing()

# Find and include aocommon git submodule
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Include aocommon
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")

find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIRS})

find_package(Casacore COMPONENTS casa ms tables fits REQUIRED)
include_directories(${CASACORE_INCLUDE_DIR})

find_package(Threads REQUIRED)

find_package(IDGAPI NO_MODULE QUIET)

if(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  include_directories(${IDGAPI_INCLUDE_DIRS})
  set(HAVE_IDG TRUE)
  add_definitions(-DHAVE_IDG)
  message(STATUS "Image domain gridder API libraries found.")
else(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  set(IDGAPI_LIBRARIES)
  message(WARNING "Image domain gridder API libraries NOT found. IDG prediction will not be available.")
endif(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})

find_package(Python3 COMPONENTS Interpreter Development)
message(STATUS "Using python version ${Python3_VERSION}")

include_directories(${Python3_INCLUDE_DIRS})

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
set(BOOSTPY ${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
find_package(Boost COMPONENTS date_time filesystem python${BOOSTPY} numpy${BOOSTPY} program_options system unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_path(EVERYBEAM_INCLUDE_DIR NAMES EveryBeam/load.h)
if(NOT EVERYBEAM_INCLUDE_DIR)
  message(FATAL_ERROR "EveryBeam include file (EveryBeam/load.h) not found: please install https://git.astron.nl/RD/EveryBeam")
endif()
find_library(EVERYBEAM_LIB everybeam)
if(EVERYBEAM_INCLUDE_DIR AND EVERYBEAM_LIB)
  include_directories(${EVERYBEAM_INCLUDE_DIR})
  message(STATUS "EveryBeam library found.")
else()
  message(FATAL_ERROR "EveryBeam library not found: please install https://git.astron.nl/RD/EveryBeam")
endif(EVERYBEAM_INCLUDE_DIR AND EVERYBEAM_LIB)

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
  DPPP/Apply.cc
  DPPP/ApplyCal.cc
  DPPP/Averager.cc
  DPPP/BaselineSelection.cc
  DPPP/BDAAverager.cc
  DPPP/BDABuffer.cc
  DPPP/BDAIntervalBuffer.cc
  DPPP/Counter.cc
  DPPP/Demixer.cc
  DPPP/DPBuffer.cc
  DPPP/DPInfo.cc
  DPPP/DPInput.cc
  DPPP/DPLogger.cc
  DPPP/DPRun.cc
  DPPP/DPStep.cc
  DPPP/DummyStep.cc
  DPPP/EstimateMixed.cc
  DPPP/EstimateNew.cc
  DPPP/Filter.cc
  DPPP/FlagCounter.cc
  DPPP/GainCal.cc
  DPPP/GainCalAlgorithm.cc
  DPPP/GaussianSource.cc
  DPPP/GridInterpolate.cc
  DPPP/H5Parm.cc
  DPPP/H5ParmPredict.cc
  DPPP/Interpolate.cc
  DPPP/MedFlagger.cc
  DPPP/ModelComponent.cc
  DPPP/ModelComponentVisitor.cc
  DPPP/MSReader.cc
  DPPP/MSUpdater.cc
  DPPP/MSWriter.cc
  DPPP/MultiMSReader.cc
  DPPP/OneApplyCal.cc
  DPPP/Patch.cc
  DPPP/PhaseFitter.cc
  DPPP/PhaseShift.cc
  DPPP/PointSource.cc
  DPPP/Position.cc
  DPPP/Predict.cc
  DPPP/PreFlagger.cc
  DPPP/ProgressMeter.cc
  DPPP/ScaleData.cc
  DPPP/SetBeam.cc
  DPPP/Simulate.cc
  DPPP/Simulator.cc
  DPPP/SolTab.cc
  DPPP/SourceDBUtil.cc
  DPPP/Split.cc
  DPPP/StationAdder.cc
  DPPP/Stokes.cc
  DPPP/SubtractMixed.cc
  DPPP/SubtractNew.cc
  DPPP/Upsample.cc
  DPPP/UVWCalculator.cc
  DPPP/UVWFlagger.cc
  DPPP/ApplyBeam.cc
  DPPP/DemixerNew.cc
  DPPP/DemixInfo.cc
  DPPP/DemixWorker.cc
)
set(DPPP_OBJECT $<TARGET_OBJECTS:DPPP_OBJ>)

find_package(AOFlagger 3 REQUIRED)
include_directories(${AOFLAGGER_INCLUDE_DIR})

add_library(AOFlaggerStep_OBJ OBJECT
  AOFlaggerStep/AOFlaggerStep.cc
  AOFlaggerStep/Register.cc
)
set(AOFLAGGERSTEP_OBJECT $<TARGET_OBJECTS:AOFlaggerStep_OBJ>)

set(EXTRA_LIBRARIES
  ${AOFLAGGER_LIB}
  ${CASACORE_LIBRARIES}
  ${Boost_LIBRARIES}
  ${EXTRA_LIBRARIES}
  ${Python3_LIBRARIES}
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

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  add_library(dppp_testdyndppp SHARED
    TestDynDPPP/TestDynStep.cc ${DPPP_OBJECT} ${PARMDB_OBJECT} ${BLOB_OBJECT} ${DDECAL_OBJECT} ${COMMON_OBJECT} ${AOFLAGGERSTEP_OBJECT}
  )
  set_target_properties(dppp_testdyndppp PROPERTIES VERSION ${DPPP_GLOBAL_SOVERSION})
endif()

# I'm avoiding cmake's test commands, because they build executable tests by default
# which makes the compilation process more expensive.
set(TEST_FILENAMES
  DPPP/test/runtests.cc
  DPPP/test/unit/mock/MockStep.cc
  DPPP/test/unit/tApplyCal.cc
  DPPP/test/unit/tApplyCalH5.cc
  DPPP/test/unit/tAverager.cc
  DPPP/test/unit/tBaselineSelection.cc
  DPPP/test/unit/tBDABuffer.cc
  DPPP/test/unit/tBDAIntervalBuffer.cc
  # DPPP/test/unit/tDemixer.cc # Parset is no longer valid in this test
  DPPP/test/unit/tFilter.cc
  DPPP/test/unit/tGridInterpolate.cc
  # DPPP/test/unit/tH5Parm.cc # Failes because sol001 cannot be found
  DPPP/test/unit/tInterpolate.cc
  DPPP/test/unit/tMedFlagger.cc
  DPPP/test/unit/tMirror.cc
  # DPPP/test/unit/tNDPPP.cc # Integration test, not yet working
  DPPP/test/unit/tPreFlagger.cc
  DPPP/test/unit/tPhaseShift.cc
  DPPP/test/unit/tUpsample.cc
  DPPP/test/unit/tUVWFlagger.cc
  DPPP/test/unit/tPSet.cc
  DPPP/test/unit/tScaleData.cc
  DPPP/test/unit/tScaleDataBDA.cc
  DPPP/test/unit/tStationAdder.cc
  DDECal/test/unit/tRotationConstraint.cc
)

# Boost 1.59 introduced BOOST_TEST. The tests below use this feature.
if(Boost_VERSION_STRING VERSION_GREATER_EQUAL "1.59")
  list(APPEND TEST_FILENAMES
    DPPP/test/unit/tBDAAverager.cc
    DPPP/test/unit/tDPInfo.cc
  )
else()
  message(WARNING "Boost < 1.59 detected. Some unit tests are disabled.")
endif()

# Add boost dynamic link flag for all test files.
# https://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/usage_variants.html
# Without this flag, linking is incorrect and boost performs duplicate delete()
# calls after running all tests, in the cleanup phase.
set_source_files_properties(
  ${TEST_FILENAMES} PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK"
)

add_test(
  NAME extract_applycal
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tApplyCal_tmp.parmdb.tgz
)
set_tests_properties(
  extract_applycal
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP applycal
)

add_test(
  NAME extract_NDPPP
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tNDPPP_tmp.MS.tgz
)
set_tests_properties(
  extract_NDPPP
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP ndppp
)

add_executable(
  unittests EXCLUDE_FROM_ALL
  ${TEST_FILENAMES} ${OS_SPECIFIC_TESTS}
  ${DPPP_OBJECT} ${COMMON_OBJECT} ${BLOB_OBJECT} ${DDECAL_OBJECT}
  ${IDGPREDICT_OBJECT} ${PARMDB_OBJECT} ${AOFLAGGERSTEP_OBJECT}
)
target_link_libraries(unittests ${EXTRA_LIBRARIES} ${EVERYBEAM_LIB})

add_test(buildunittests ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unittests)
set_tests_properties(buildunittests PROPERTIES FIXTURES_SETUP unittests)

add_test(
  NAME unittests
  COMMAND unittests -f JUNIT -k unittests.xml --catch_system_error=yes
)
set_tests_properties(
  unittests PROPERTIES LABELS unit
  FIXTURES_REQUIRED "unittests;applycal"
)

# I'll leave this here as an example how separate test executables could be added
#add_executable(tMirror EXCLUDE_FROM_ALL
#  DPPP/test/tMirror.cc
#)
#target_link_libraries(tMirror ${CASACORE_LIBRARIES})

install (TARGETS makesourcedb showsourcedb DESTINATION bin)

include(docs/CMakeLists.txt)
