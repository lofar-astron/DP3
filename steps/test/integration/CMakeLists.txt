# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

configure_file(testInit.sh.in testInit.sh)

# The 'source' symbolic link simplifies running the tests manually inside
# ${CMAKE_CURRENT_BINARY_DIR}: It allows using 'source/tApplyBeam.sh' instead of
# '../../../../DP3/steps/test/integration/tApplyBeam.sh. (Using 'RESULT', fatal
# errors won't occur on systems without symlink support.)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14")
  file(
    CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/source
    RESULT DUMMY_RESULT
    SYMBOLIC)
endif()

set(INTEGRATION_TESTS
    tApplyBeam
    tApplyCal2
    tBdaPredict
    tDemix
    tGainCal
    tGainCalH5Parm
    tMultiApplyCal
    tPredict
    tReadOnly
    tSplit)

foreach(TEST ${INTEGRATION_TESTS})
  # Use ${CMAKE_CURRENT_SOURCE_DIR} instead of 'source' since not all systems
  # support symbolic links.
  add_test(${TEST} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.sh)
  set_tests_properties(${TEST} PROPERTIES LABELS "integration")
endforeach()
