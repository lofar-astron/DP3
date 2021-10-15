# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../CMake)

project(DP3_DOC NONE)

find_package(Doxygen)

if(NOT DP3_VERSION)
  # Building docs standalone, get version from top level CMakelists
  file(READ "../CMakeLists.txt" TOPLEVEL_CMAKELISTS)
  set(DP3_VERSION_REGEX
      ".*set\\(DP3_VERSION ([0-9]+)\\.([0-9]+)\\.([0-9]+).*\\)")
  if("${TOPLEVEL_CMAKELISTS}" MATCHES ${DP3_VERSION_REGEX})
    set(DP3_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set(DP3_VERSION_MINOR "${CMAKE_MATCH_2}")
    set(DP3_VERSION_PATCH "${CMAKE_MATCH_3}")
  else()
    message(
      FATAL_ERROR "Failed to parse DP3_VERSION from top level CMakeLists.txt")
  endif()
endif()

if(DOXYGEN_FOUND)
  configure_file(doxygen/Doxyfile.in Doxyfile @ONLY)
  add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif(DOXYGEN_FOUND)

find_package(Sphinx)
find_program(
  JSONSCHEMA2RST_EXECUTABLE
  NAMES jsonschema2rst
  DOC "Path to jsonschema2rst executable")

if(SPHINX_FOUND AND JSONSCHEMA2RST_EXECUTABLE)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR})
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/)
  configure_file(index.rst index.rst COPYONLY)
  configure_file(conf.py.in conf.py)
  add_custom_target(
    userdoc
    COMMAND ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating user documentation with Sphinx")
else()
  message(
    "Sphinx (sphinx-build) and jsonschema2rst are necessary for building user docs"
  )
endif(SPHINX_FOUND AND JSONSCHEMA2RST_EXECUTABLE)
