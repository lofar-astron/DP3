# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

pybind11_add_module(pydppp MODULE pydppp.cc)
target_link_libraries(pydppp PUBLIC ${CASACORE_LIBRARIES})

pybind11_add_module(parameterset MODULE parameterset.cc)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--unresolved-symbols=ignore-all")

install(TARGETS pydppp parameterset DESTINATION "${PYTHON_INSTALL_DIR}/dppp")
install(FILES __init__.py dummystep.py DESTINATION "${PYTHON_INSTALL_DIR}/dppp")