# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later


pybind11_add_module(pydppp MODULE pydppp.cc)

pybind11_add_module(parameterset MODULE parameterset.cc)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--unresolved-symbols=ignore-all")

add_library(dppp_pythondppp SHARED PyDPStep.cc)

target_link_libraries(dppp_pythondppp pybind11::embed)

install(TARGETS dppp_pythondppp)
install(TARGETS pydppp parameterset DESTINATION "${PYTHON_INSTALL_DIR}/dppp")
install(FILES __init__.py dummystep.py DESTINATION "${PYTHON_INSTALL_DIR}/dppp")


#  add_subdirectory(test)

# Commented out code below should move to test subdirectory

#add_executable(testpydpstep testmydpstep.cc
#  ${DPPP_OBJECT}
#  ${PARMDB_OBJECT}
#  ${BLOB_OBJECT}
#  ${COMMON_OBJECT}
#  ${AOFLAGGERSTEP_OBJECT}
#  ${DDECAL_OBJECT}
#  ${IDGPREDICT_OBJECT}
#)
#target_link_libraries(testpydpstep PRIVATE PythonStep)
#set_target_properties(testpydpstep PROPERTIES ENABLE_EXPORTS ON)

#target_link_libraries(testpydpstep PRIVATE
#  ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}
#  ${CASACORE_LIBRARIES}
#  ${EVERYBEAM_LIB}
#  ${CFITSIO_LIBRARY}
#  ${EXTRA_LIBRARIES})