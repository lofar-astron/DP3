# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Ensure that we use /usr as install prefix when packaging.
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")

# We use git to determine the version string from the latest tag.
find_package(Git)
if(Git_FOUND)
  execute_process(
    COMMAND bash -c "${GIT_EXECUTABLE} describe --tags --dirty"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Strip leading non-numeric characters, like "v"
  string(REGEX REPLACE "^[^[0-9]]*" "" VERSION_STRING ${VERSION_STRING})
endif()
if(NOT VERSION_STRING)
  set(VERSION_STRING "0.0")
  message(WARNING "[CPack] Could not find tag to create version string, "
                  "using ${VERSION_STRING} instead.")
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/DetermineTargetCPU.cmake)

# Set package name: encode target CPU, if defined, in package name.
if(DEFINED IDENTIFIED_TARGET_CPU)
  set(CPACK_PACKAGE_NAME "dp3-${IDENTIFIED_TARGET_CPU}")
else()
  set(CPACK_PACKAGE_NAME "dp3")
endif()

set(CPACK_PACKAGE_VENDOR "ASTRON")
set(CPACK_PACKAGE_VERSION "${VERSION_STRING}")

set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LOFAR preprocessing software")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://git.astron.nl/RD/DP3")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")

set(CPACK_GENERATOR "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "aoflagger (>= 3.0.1),\
 everybeam (>= 0.2), everybeam (<< 0.3),\
 idg-api (>= 0.8)")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "deb-packages@astron.nl")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# Determine list of conflicting package names. Packages not built for current
# target CPU always conflict.
if(DEFINED IDENTIFIED_TARGET_CPU)
  set(_conflicts dp3)
else()
  set(_conflicts)
endif()
foreach(_cpu ${KNOWN_TARGET_CPUS})
  if(NOT "${_cpu}" STREQUAL "${IDENTIFIED_TARGET_CPU}")
    if("${_conflicts}" STREQUAL "")
      set(_conflicts "dp3-${_cpu}")
    else()
      set(_conflicts "${_conflicts}, dp3-${_cpu}")
    endif()
  endif()
endforeach()
set(CPACK_DEBIAN_PACKAGE_CONFLICTS "${_conflicts}")

include(CPack)

message(STATUS "Package name: ${CPACK_PACKAGE_NAME}")
message(STATUS "Package version: ${CPACK_PACKAGE_VERSION}")
