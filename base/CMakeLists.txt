# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(DP3 Main.cc ${DP3_OBJECTS})
target_link_libraries(DP3 ${DP3_LIBRARIES})
add_dependencies(DP3 schaapcommon)

set_target_properties(DP3 PROPERTIES ENABLE_EXPORTS ON)
set_target_properties(DP3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                     ${CMAKE_BINARY_DIR})

install(TARGETS DP3 DESTINATION bin)

# Data files
install(FILES share/LOFAR-LBA-default.rfis share/LOFAR-HBA-default.rfis
        DESTINATION share/rfistrategies)
