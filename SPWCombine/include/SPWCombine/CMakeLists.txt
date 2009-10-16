# $Id$

# List of header files that will be installed.
set(inst_HEADERS
  CombinerProcessControl.h
  SPWCombine.h)

# Add definition of include path suffix definition, because the header files
# reside in a non-standard location.
set(${PACKAGE_NAME}_INCLUDE_PATH_SUFFIX SPWCombine
  CACHE INTERNAL "Include path suffix for package \"${PACKAGE_NAME}\"")

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/SPWCombine)

# Install header files.
install(FILES ${inst_HEADERS} DESTINATION include/${PACKAGE_NAME})
