configure_file(testInit.sh.in testInit.sh)

# The 'source' symbolic link simplifies running the tests manually inside
# ${CMAKE_CURRENT_BINARY_DIR}: It allows using 'source/tApplyBeam.sh' instead of
# '../../../../DP3/DPPP/test/integration/tApplyBeam.sh.
# (Using 'RESULT', fatal errors won't occur on systems without symlink support.)
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/source
     RESULT DUMMY_RESULT SYMBOLIC)

# Use ${CMAKE_CURRENT_SOURCE_DIR} instead of 'source' since not all systems
# support symbolic links.
add_test(tApplyBeam ${CMAKE_CURRENT_SOURCE_DIR}/tApplyBeam.sh)
add_test(tGainCal ${CMAKE_CURRENT_SOURCE_DIR}/tGainCal.sh)
add_test(tGainCalH5Parm ${CMAKE_CURRENT_SOURCE_DIR}/tGainCalH5Parm.sh)
add_test(tPredict ${CMAKE_CURRENT_SOURCE_DIR}/tPredict.sh)
