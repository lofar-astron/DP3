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

set(INTEGRATION_TESTS tBdaPredict tDemix tGainCal)

foreach(TEST ${INTEGRATION_TESTS})
  # Use ${CMAKE_CURRENT_SOURCE_DIR} instead of 'source' since not all systems
  # support symbolic links.
  add_test(${TEST} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.sh)
  set_tests_properties(${TEST} PROPERTIES LABELS "integration")
endforeach()

# Example pytest integration tests
# Notes:
# - python / pytest testing script will replace the bash integration tests
# - the collection of python-based integration tests can be invoked with ctest -L pyintegration
# - see tPredict.py docstring how the tPredict.py tests can be invoked stand-alone
set(PYINTEGRATION_TESTS
    tApplyBeam
    tApplyCal2
    tBdaExpander
    tGainCalH5Parm
    tMultiApplyCal
    tPredict
    tReadOnly
    tSplit)

foreach(TEST ${PYINTEGRATION_TESTS})
  add_test(
    NAME ${TEST}
    COMMAND
      pytest --exitfirst -v
      --junitxml=${CMAKE_BINARY_DIR}/integration_${TEST}.xml
      "${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.py")
  set_tests_properties(${TEST} PROPERTIES LABELS "pyintegration")
endforeach()
