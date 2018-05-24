# $Id: CMakeLists.txt 39071 2018-02-08 15:50:02Z dijkema $

# List of header files that will be installed.
set(inst_HEADERS
  buffered_lane.h Interpolate.h lane.h
  )

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/${PACKAGE_NAME})

# Install header files.
install(FILES ${inst_HEADERS} DESTINATION include/${PACKAGE_NAME})
