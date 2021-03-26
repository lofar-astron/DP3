# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

pybind11_add_module(pydp3 MODULE pydp3.cc)
target_link_libraries(pydp3 PUBLIC ${CASACORE_LIBRARIES})

pybind11_add_module(parameterset MODULE parameterset.cc)

install(TARGETS pydp3 parameterset DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
install(FILES __init__.py DESTINATION "${PYTHON_INSTALL_DIR}/dp3")
