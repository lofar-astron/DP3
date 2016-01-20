# $Id: CMakeLists.txt 30990 2015-02-12 12:27:47Z diepen $

# List of header files that will be installed.
set(inst_HEADERS AOFlaggerStep.h SlidingFlagger.h Register.h)

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/${PACKAGE_NAME})

# Install header files.
#install(FILES ${inst_HEADERS} DESTINATION include/${PACKAGE_NAME})
