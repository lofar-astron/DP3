# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Create a pseudo install directory 'dp3' in the build directory that
# mimics the ${PYTHON_INSTALL_DIR}$ such that 'import dp3' can be used
# in the build environment
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/dp3")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
               "${CMAKE_CURRENT_BINARY_DIR}/dp3/__init__.py" COPYONLY)

pybind11_add_module(pydp3 MODULE pydp3.cc pydpinfo.cc)
target_link_libraries(pydp3 PRIVATE LIBDP3)
# Output the library/module in the pseudo install directory 'dp3'
set_target_properties(pydp3 PROPERTIES LIBRARY_OUTPUT_DIRECTORY dp3)

pybind11_add_module(parameterset MODULE parameterset.cc)
target_link_libraries(parameterset PRIVATE LIBDP3)
# Output the library/module in the pseudo install directory 'dp3'
set_target_properties(parameterset PROPERTIES LIBRARY_OUTPUT_DIRECTORY dp3)

install(TARGETS parameterset pydp3 DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
install(FILES __init__.py DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
