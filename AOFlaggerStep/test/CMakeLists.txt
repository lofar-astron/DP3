# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Do not use lofar_add_test, because libaoflaggerstep should not be linked in.
add_executable(tAOFlaggerStep tAOFlaggerStep.cc)

target_link_libraries (tAOFlaggerStep ${_libs}  ${LOFAR_EXTRA_LIBRARIES})

add_test (tAOFlaggerStep ${CMAKE_SOURCE_DIR}/DP3/DPPP_AOFlag/test/tAOFlaggerStep.sh)
