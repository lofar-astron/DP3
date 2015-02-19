# $Id: CMakeLists.txt 26355 2013-09-04 11:51:07Z dijkema $

include(LofarCTest)

# Do not use lofar_add_test, because libtestdyndppp should not be linked in.
add_executable (tTestDynStep tTestDynStep.cc)
target_link_libraries (tTestDynStep ${_libs}  ${LOFAR_EXTRA_LIBRARIES})
add_test (tTestDynStep ${CMAKE_SOURCE_DIR}/CEP/DP3/TestDynDPPP/test/tTestDynStep.sh)
