# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# The .cc tests are added in the top level CMakeLists.txt file.

if("${CMAKE_CXX_FLAGS}" MATCHES "-fsanitize=address")
  message(STATUS "The Python tests do not work with address sanitizer enabled.")
else()
  add_python_tests(
    tPyDpBuffer
    tPyDpInfo
    tPyFields
    tPyFitters
    tPyParameterSet
    tPyStep
    tQueueOutput)
endif()
