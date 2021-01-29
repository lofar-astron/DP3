# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# FindHDF5 uses NATIVE_COMMAND in separate_arguments, which requires
# CMake 3.9.
cmake_minimum_required(VERSION 3.9)

# CMake >= 3.19.1 gives a warning when this policy is not 'NEW'.
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.19.1")
  cmake_policy(SET CMP0074 NEW)
endif()

# Set version number and project name.
set(DP3_VERSION 5.0.0)
if(DP3_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)")
  set(DP3_VERSION_MAJOR "${CMAKE_MATCH_1}")
  set(DP3_VERSION_MINOR "${CMAKE_MATCH_2}")
  set(DP3_VERSION_PATCH "${CMAKE_MATCH_3}")
else()
  message(FATAL_ERROR "Failed to parse DP3_VERSION='${DP3_VERSION}'")
endif()

project(DP3 VERSION ${DP3_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

add_compile_options(-Wall -O3)
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")

if(NOT PORTABLE)
  add_compile_options(-march=native)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)


# === Load external packages ===

# DDECal dependencies
find_package(Armadillo)
if(${ARMADILLO_FOUND})
  add_definitions(-DHAVE_ARMADILLO)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  set(DDE_ARMADILLO_FILES
    DDECal/KLFitter.cc
    DDECal/PiercePoint.cc
    DDECal/ScreenConstraint.cc)
else()
  message(WARNING "Armadillo was not found, not including screenfitter inside DDECal")
  set(DDE_ARMADILLO_FILES)
endif()

# Casacore depends on HDF5 -> First load HDF5.
find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIRS})

find_package(Casacore COMPONENTS casa ms tables fits REQUIRED)
include_directories(${CASACORE_INCLUDE_DIR})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12")
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  set(BOOSTPY ${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
else() # Use old, deprecated means of detecting python.
  find_package(PythonInterp 3 REQUIRED)
  find_package(PythonLibs 3 REQUIRED)
  set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
  set(Python3_VERSION ${PYTHON_VERSION_STRING})
  set(Python3_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
  set(Python3_LIBRARIES ${PYTHON_LIBRARIES})
endif()
message(STATUS "Using python version ${Python3_VERSION}")
include_directories(${Python3_INCLUDE_DIRS})

#Boost depends on python, since we use ${BOOSTPY} here.
#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time filesystem program_options system unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_package(Threads REQUIRED)


# === Load astron packages ===

find_package(AOFlagger 3.0.1 REQUIRED)
include_directories(${AOFLAGGER_INCLUDE_DIR})

# We could use find_package(EveryBeam 0.1.1 REQUIRED),
# however conditions below make it somewhat more explicit
find_package(EveryBeam NO_MODULE)
if(${EVERYBEAM_FOUND})
  if(${EVERYBEAM_VERSION} VERSION_LESS "0.1.1" OR ${EVERYBEAM_VERSION} VERSION_GREATER_EQUAL "0.2.0")
    message(FATAL_ERROR "DP3 needs EveryBeam version 0.1.x - with x > 0 - but found version ${EveryBeam_VERSION}")
  endif()
include_directories(${EVERYBEAM_INCLUDE_DIR})
else(${EVERYBEAM_FOUND})
  message(FATAL_ERROR "DP3 requires EveryBeam, but EveryBeam was not found. "
    "Please install https://git.astron.nl/RD/EveryBeam and make sure that "
    "EveryBeam is added to the CMAKE_PREFIX_PATH")
endif(${EVERYBEAM_FOUND})

find_package(IDGAPI NO_MODULE QUIET)
if(IDGAPI_FOUND)
  # Throw error if IDG version < 0.8 or version not provided at all
  if((IDGAPI_VERSION VERSION_LESS "0.8") OR (NOT DEFINED IDGAPI_VERSION))
    message(FATAL_ERROR "IDGAPI was found, but DP3 requires IDGAPI to have version >= 0.8. "
    "Please compile IDG repository at a version >= 0.8")
  endif()
endif()
if(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  include_directories(${IDGAPI_INCLUDE_DIRS})
  set(HAVE_IDG TRUE)
  add_definitions(-DHAVE_IDG)
  message(STATUS "Image domain gridder API libraries found.")
else(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  set(IDGAPI_LIBRARIES "")
  message(WARNING "Image domain gridder API libraries NOT found. IDG prediction will not be available.")
endif(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)


# === Load internal submodule packages. ===

# Update submodules as needed
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Update submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Syncing submodules")
    # Account for potential changes in git repo URL's
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule sync --recursive
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "Syncing submodules did not succeed"
              "command '${GIT_EXECUTABLE} submodule sync --recursive' failed with exit code ${GIT_SUBMOD_RESULT}")
    endif()
    message(STATUS "Updating submodules")
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

# Include pybind11
set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(${pybind11_INCLUDE_DIR})

# Include schaapcommon
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")

# Add cmake information to headers
configure_file(DPPP/Version.h.in DPPP/Version.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/DPPP)

# Object libraries (which are just file groups, not actual libraries) are made in order
# to avoid compiling everything multiple times for multiple executables.
add_library(AOFlaggerStep OBJECT
  AOFlaggerStep/AOFlaggerStep.cc
  AOFlaggerStep/Register.cc
)

add_library(Blob OBJECT
  Blob/BlobAipsIO.cc
  Blob/BlobArray.cc
  Blob/BlobHeader.cc
  Blob/BlobIBufStream.cc
  Blob/BlobIStream.cc
  Blob/BlobOBufStream.cc
  Blob/BlobOStream.cc
)

add_library(Common OBJECT
  Common/BaselineSelect.cc
  Common/ClusterDesc.cc
  Common/DataConvert.cc
  Common/Memory.cc
  Common/NodeDesc.cc
  Common/ParameterHandler.cc
  Common/ParameterRecord.cc
  Common/ParameterSet.cc
  Common/ParameterSetImpl.cc
  Common/ParameterValue.cc
  Common/PrettyUnits.cc
  Common/ProximityClustering.cc
  Common/StringUtil.cc
  Common/Timer.cc
  Common/TypeNames.cc
  Common/VdsDesc.cc
  Common/VdsMaker.cc
  Common/VdsPartDesc.cc
)

add_library(DDECal OBJECT
  DDECal/Constraint.cc
  DDECal/DDECal.cc
  DDECal/DiagonalSolver.cc
  DDECal/FullJonesSolver.cc
  DDECal/IterativeDiagonalSolver.cc
  DDECal/RotationAndDiagonalConstraint.cc
  DDECal/RotationConstraint.cc
  DDECal/ScalarSolver.cc
  DDECal/SmoothnessConstraint.cc
  DDECal/SolverBase.cc
  DDECal/Stopwatch.cc
  DDECal/TECConstraint.cc
  ${DDE_ARMADILLO_FILES})

add_library(DPPP_OBJ OBJECT
  DPPP/Apply.cc
  DPPP/ApplyCal.cc
  DPPP/Averager.cc
  DPPP/BaselineSelection.cc
  DPPP/BDAAverager.cc
  DPPP/BDABuffer.cc
  DPPP/BDAIntervalBuffer.cc
  DPPP/ColumnReader.cc
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
  DPPP/H5ParmPredict.cc
  DPPP/Interpolate.cc
  DPPP/MedFlagger.cc
  DPPP/ModelComponent.cc
  DPPP/ModelComponentVisitor.cc
  DPPP/MS.cc
  DPPP/MSBDAReader.cc
  DPPP/MSBDAWriter.cc
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
  DPPP/SolutionInterval.cc
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

add_library(IDGPredict OBJECT
  IDGPredict/IDGPredict.cc
)

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

add_library(PythonDPPP OBJECT
  PythonDPPP/PyDPStep.cc
)

# dppp_testdyndppp requires using position independent code.
set_property(TARGET AOFlaggerStep Blob Common DDECal DPPP_OBJ IDGPredict ParmDB PythonDPPP
             PROPERTY POSITION_INDEPENDENT_CODE ON)

set(DPPP_OBJECTS
  $<TARGET_OBJECTS:AOFlaggerStep>
  $<TARGET_OBJECTS:Blob>
  $<TARGET_OBJECTS:Common>
  $<TARGET_OBJECTS:DDECal>
  $<TARGET_OBJECTS:DPPP_OBJ>
  $<TARGET_OBJECTS:IDGPredict>
  $<TARGET_OBJECTS:ParmDB>
  $<TARGET_OBJECTS:PythonDPPP>)

set(DPPP_LIBRARIES
  ${AOFLAGGER_LIB}
  ${ARMADILLO_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CASACORE_LIBRARIES}
  ${CFITSIO_LIBRARY}
  ${EVERYBEAM_LIB}
  ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}
  ${IDGAPI_LIBRARIES}
  ${Python3_LIBRARIES}
  schaapcommon
  Threads::Threads
  pybind11::embed)

add_subdirectory(DPPP)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  add_subdirectory(PythonDPPP)
endif()
# SPW Combine has dependencies to LOFAR
#add_subdirectory(SPWCombine)

set(SOURCEDB_OBJECTS
  $<TARGET_OBJECTS:Blob>
  $<TARGET_OBJECTS:Common>
  $<TARGET_OBJECTS:ParmDB>)

set(SOURCEDB_LIBRARIES
  ${CASACORE_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY})

add_executable(makesourcedb ParmDB/makesourcedb.cc ${SOURCEDB_OBJECTS})
target_link_libraries(makesourcedb ${SOURCEDB_LIBRARIES})

add_executable(showsourcedb ParmDB/showsourcedb.cc ${SOURCEDB_OBJECTS})
target_link_libraries(showsourcedb ${SOURCEDB_LIBRARIES})

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(OS_SPECIFIC_TESTS) # No specific tests for Apple
else()
  # These run only on Linux
  set(OS_SPECIFIC_TESTS
#    TestDynDPPP/test/tDynStep.cc # This test still fails
  )
  add_library(dppp_testdyndppp SHARED
    TestDynDPPP/TestDynStep.cc ${DPPP_OBJECTS}
  )
  target_link_libraries(dppp_testdyndppp ${DPPP_LIBRARIES})
  add_dependencies(dppp_testdyndppp schaapcommon)
endif()

# I'm avoiding cmake's test commands, because they build executable tests by default
# which makes the compilation process more expensive.
set(TEST_FILENAMES
  Common/test/unit/tProximityClustering.cc
  DPPP/test/runtests.cc
  DPPP/test/unit/mock/MockInput.cc
  DPPP/test/unit/mock/MockStep.cc
  DPPP/test/unit/fixtures/fDirectory.cc
  DPPP/test/unit/tApplyCal.cc
  DPPP/test/unit/tApplyCalH5.cc
  DPPP/test/unit/tAverager.cc
  DPPP/test/unit/tBaselineSelection.cc
  DPPP/test/unit/tBDABuffer.cc
  DPPP/test/unit/tBDAIntervalBuffer.cc
  # DPPP/test/unit/tDemixer.cc # Parset is no longer valid in this test
  DPPP/test/unit/tFilter.cc
  DPPP/test/unit/tInterpolate.cc
  DPPP/test/unit/tMedFlagger.cc
  DPPP/test/unit/tDPPP.cc
  DPPP/test/unit/tMirror.cc
  DPPP/test/unit/tMSReader.cc
  DPPP/test/unit/tPreFlagger.cc
  DPPP/test/unit/tPhaseShift.cc
  DPPP/test/unit/tUpsample.cc
  DPPP/test/unit/tUVWFlagger.cc
  DPPP/test/unit/tPyDPStep.cc
  DPPP/test/unit/tPSet.cc
  DPPP/test/unit/tScaleData.cc
  DPPP/test/unit/tSolutionInterval.cc
  DPPP/test/unit/tSourceDBUtil.cc
  DPPP/test/unit/tStationAdder.cc
  DDECal/test/unit/tRotationConstraint.cc
)
set(TEST_SOLVER_FILENAMES
  DPPP/test/runtests.cc
  DDECal/test/unit/tSolvers.cc
)

# Boost 1.59 introduced BOOST_TEST. The tests below use this feature.
if(Boost_VERSION_STRING VERSION_GREATER_EQUAL "1.59")
  list(APPEND TEST_FILENAMES
    Common/test/unit/tMemory.cc
    Common/test/unit/tStringUtil.cc
    DPPP/test/unit/tBDAAverager.cc
    DPPP/test/unit/tDPInfo.cc
    DPPP/test/unit/tDPInput.cc
    DPPP/test/unit/tMSBDAReader.cc
    DPPP/test/unit/tMSBDAWriter.cc
    DPPP/test/unit/tScaleDataBDA.cc)
  if(HAVE_IDG)
    list(APPEND TEST_FILENAMES IDGPredict/test/unit/tIDGPredict.cc)
  endif()
else()
  message(WARNING "Boost < 1.59 detected. Some unit tests are disabled.")
endif()

# Add boost dynamic link flag for all test files.
# https://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/usage_variants.html
# Without this flag, linking is incorrect and boost performs duplicate delete()
# calls after running all tests, in the cleanup phase.
set_source_files_properties(
  ${TEST_FILENAMES} ${TEST_SOLVER_FILENAMES} PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK"
)

add_test(
  NAME extract_applycal
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tApplyCal_tmp.parmdb.tgz
)
set_tests_properties(
  extract_applycal
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_test(
  NAME extract_NDPPP
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tNDPPP.in_MS.tgz
)
set_tests_properties(
  extract_NDPPP
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_test(
  NAME extract_NDPPP_BDA
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tNDPPP_bda.in_MS.tgz
)
set_tests_properties(
  extract_NDPPP_BDA
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_test(
  NAME extract_NDPPP_generic
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/integration/tNDPPP-generic.MS.tgz
)
set_tests_properties(
  extract_NDPPP_generic
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_test(
  NAME extract_OSKAR
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/DPPP/test/resources/tOSKAR.in_MS.tgz
)
set_tests_properties(
  extract_OSKAR
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_test(
  NAME extract_IDG
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/IDGPredict/test/resources/sources.tgz
)
set_tests_properties(
  extract_IDG
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} FIXTURES_SETUP extract_resources
)

add_executable(
  unittests EXCLUDE_FROM_ALL
  ${TEST_FILENAMES} ${OS_SPECIFIC_TESTS} ${DPPP_OBJECTS}
)
set_target_properties(unittests PROPERTIES ENABLE_EXPORTS ON)
target_link_libraries(unittests ${DPPP_LIBRARIES})
add_dependencies(unittests schaapcommon)

add_executable(
  unittest_solver EXCLUDE_FROM_ALL
  ${TEST_SOLVER_FILENAMES} ${OS_SPECIFIC_TESTS} ${DPPP_OBJECTS}
)
target_link_libraries(unittest_solver ${DPPP_LIBRARIES})

add_test(buildunittests ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unittests unittest_solver)
set_tests_properties(buildunittests PROPERTIES FIXTURES_SETUP unittests)
add_subdirectory(DDECal/test/integration)

# unittests.sh adjusts the PYTHONPATH to make tPyDPStep working.
configure_file(unittests.sh.in unittests.sh)
add_test(
  NAME unittests
  COMMAND unittests.sh -f JUNIT -k unittests.xml --catch_system_error=yes
)
set_tests_properties(
  unittests PROPERTIES LABELS unit
  FIXTURES_REQUIRED "unittests;extract_resources"
)

add_test(
  NAME unittest_solver
  COMMAND unittest_solver -f JUNIT -k unittest_solver.xml --catch_system_error=yes
)
set_tests_properties(
  unittest_solver PROPERTIES LABELS unit
  FIXTURES_REQUIRED "unittests"
)

install(TARGETS makesourcedb showsourcedb DESTINATION bin)

include(docs/CMakeLists.txt)
