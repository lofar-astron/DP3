# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(PythonInstall)

add_executable(DPPP DPPP.cc ${DPPP_OBJECTS})
target_link_libraries(DPPP ${DPPP_LIBRARIES})
add_dependencies(DPPP schaapcommon)

set_target_properties(DPPP PROPERTIES ENABLE_EXPORTS ON)

install(TARGETS DPPP DESTINATION bin)

# Data files
install(
  FILES share/LOFAR-LBA-default.rfis share/LOFAR-HBA-default.rfis
  DESTINATION share/rfistrategies)

# lofar_add_bin_scripts(taqlflagger)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)

add_subdirectory(test/integration)
