# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Create a pseudo install directory 'dp3' in the build directory that
# mimics the ${PYTHON_INSTALL_DIR} such that 'import dp3' can be used
# in the build environment

set(PSEUDO_INSTALL_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dp3)

file(MAKE_DIRECTORY "${PSEUDO_INSTALL_DIRECTORY}")
file(MAKE_DIRECTORY "${PSEUDO_INSTALL_DIRECTORY}/steps")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
               "${PSEUDO_INSTALL_DIRECTORY}/__init__.py" COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/steps/__init__.py"
               "${PSEUDO_INSTALL_DIRECTORY}/steps/__init__.py" COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/steps/queue.py"
               "${PSEUDO_INSTALL_DIRECTORY}/steps/queue.py" COPYONLY)

pybind11_add_module(pydp3 MODULE pydp3.cc pydpbuffer.cc pydpinfo.cc pyfields.cc)
pybind11_add_module(parameterset MODULE parameterset.cc)
pybind11_add_module(fitters MODULE pyfitters.cc)

foreach(PYLIB pydp3 parameterset fitters)
  target_link_libraries(${PYLIB} PRIVATE LIBDP3)

  # Output the library/module in the pseudo install directory.
  set_target_properties(${PYLIB} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                            ${PSEUDO_INSTALL_DIRECTORY})
endforeach()
target_link_libraries(fitters PRIVATE xtensor)

install(TARGETS parameterset pydp3 fitters
        DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
install(FILES __init__.py DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
install(FILES steps/__init__.py DESTINATION "${PYTHON_INSTALL_DIR}/dp3/steps")
install(FILES steps/queue.py DESTINATION "${PYTHON_INSTALL_DIR}/dp3/steps")
