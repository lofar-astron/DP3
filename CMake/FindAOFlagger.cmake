# This is a FindPackage script for AOFlagger.
#
# It sets the following variables:
# - AOFLAGGER_INCLUDE_DIR
# - AOFLAGGER_LIB
# - AOFLAGGER_VERSION[_MAJOR/_MINOR]
# - AOFLAGGER_FOUND
# - AOFLAGGER_ROOT_DIR (also used as hint)
# The variable AOFLAGGER_ROOT_DIR can be set to select a specific version of aoflagger.
#
# A few dependencies of AOFlagger are not handled by this script, because they are complex and
# often already included by the caller. These are: casacore, Boost, Python
# These need to be manually configured by the caller.

cmake_minimum_required(VERSION 3.0)

set(AOFLAGGER_FOUND FALSE)

# Try to find AOFlagger config file (aoflagger-config.cmake) which knows everything
# about the aoflagger installation
if(NOT AOFLAGGER_NO_AOFLAGGER_CMAKE)
  find_package(AOFLAGGER NO_MODULE QUIET)
endif()

# If no config file found, try to guess the configuration from files found
if(NOT AOFLAGGER_FOUND)
  
  find_library(AOFLAGGER_LIB aoflagger
    HINTS ${AOFLAGGER_ROOT_DIR} PATH_SUFFIXES lib)
  if(NOT AOFLAGGER_LIB)
    message(STATUS "AOFlagger library not found.")
  else()
    find_path(AOFLAGGER_INCLUDE_DIR NAMES aoflagger.h
      HINTS ${AOFLAGGER_ROOT_DIR} PATH_SUFFIXES include)
    if(NOT AOFLAGGER_INCLUDE_DIR)
      message(STATUS "AOFlagger headers not found.")
    endif()
  endif()
  
  if(AOFLAGGER_LIB AND AOFLAGGER_INCLUDE_DIR)
    execute_process(COMMAND ${AOFLAGGER_INCLUDE_DIR}/../bin/aoflagger --version OUTPUT_VARIABLE _aoflagger_result)
    if(_aoflagger_result)
      get_filename_component(AOFLAGGER_ROOT_DIR ${AOFLAGGER_INCLUDE_DIR} DIRECTORY)
      # AOFlagger --version returns something like
      # "AOFlagger 3.0-alpha", hence split space and dash:
      string(REPLACE "-" " " _aoflagger_result ${_aoflagger_result})
      string(REPLACE " " ";" _aoflagger_reslist ${_aoflagger_result})
      list(GET _aoflagger_reslist 1 AOFLAGGER_VERSION)
      string(REPLACE "." ";" _aoflagger_version_list ${AOFLAGGER_VERSION})
      list(GET _aoflagger_version_list 0 AOFLAGGER_VERSION_MAJOR)
      list(GET _aoflagger_version_list 1 AOFLAGGER_VERSION_MINOR)
      message(STATUS "Found AOFlagger version ${AOFLAGGER_VERSION_MAJOR}.${AOFLAGGER_VERSION_MINOR}: ${AOFLAGGER_LIB}")
    else()
      message(STATUS "  Executing process aoflagger failed -- AOFlagger set to NOT FOUND.")
    endif()
  endif()

  if(AOFLAGGER_LIB AND AOFLAGGER_INCLUDE_DIR AND AOFLAGGER_VERSION)
    set(_dependencies_met TRUE)
    #
    # Dependencies of AOFlagger
    #
    find_package(LibXml2)
    list(APPEND AOFLAGGER_LIB ${LibXml2_LIB})
    if(NOT LibXml2_FOUND)
      message(STATUS "AOFlagger dependency LibXml2 not found.")
      set(_dependencies_met FALSE)
    endif()
    
    find_package(PNG)
    list(APPEND AOFLAGGER_LIB ${PNG_LIBRARIES})
    if(NOT PNG_FOUND)
      message(STATUS "AOFlagger dependency PNG not found.")
      set(_dependencies_met FALSE)
    endif()
    
    find_library(FFTW3_LIB fftw3 REQUIRED)
    list(APPEND AOFLAGGER_LIB ${FFTW3_LIB})
    if(NOT FFTW3_LIB)
      message(STATUS "AOFlagger dependency FFTW3 not found.")
      set(_dependencies_met FALSE)
    endif()
    
    find_package(CFITSIO)
    if(CFITSIO_FOUND)
      list(APPEND AOFLAGGER_LIB ${CFITSIO_LIBRARY})
      list(APPEND AOFLAGGER_INCLUDE_DIR ${CFITSIO_INCLUDE_DIR})
    else()
      message(STATUS "AOFlagger dependency CFITSIO not found.")
      set(_dependencies_met FALSE)
    endif()

    find_package(PkgConfig)
    pkg_check_modules(GTKMM gtkmm-3.0>=3.0.0)
    if(GTKMM_FOUND)
      list(APPEND AOFLAGGER_LIB ${GTKMM_LIBRARIES} ${GLIBMM_LIBRARIES})
    endif(GTKMM_FOUND)

    if(_dependencies_met)
      message(STATUS "Found all of AOFlagger's dependencies.")
      SET(AOFLAGGER_FOUND TRUE)
    else()
      message(STATUS "Could not find some of AOFlagger's dependencies.")
    endif()
  else()
    message(STATUS "Could not find AOFlagger library, header or executable.")
  endif()
endif() # End manual dependency checking (this should all be set in the AOFlagger config)

# Perform version check
if(AOFLAGGER_FOUND)
  if(AOFLAGGER_VERSION VERSION_LESS AOFlagger_FIND_VERSION OR
      NOT AOFLAGGER_VERSION_MAJOR VERSION_EQUAL AOFlagger_FIND_VERSION)
    message(STATUS "Incorrect version of AOFlagger: found ${AOFLAGGER_VERSION}, need ${AOFlagger_FIND_VERSION}")
    message(STATUS "Note:  An alternative version of AOFlagger can be selected by setting the AOFLAGGER_ROOT variable.")
    set(AOFLAGGER_FOUND FALSE)
  endif()
endif()

if(NOT AOFLAGGER_FOUND)
  unset(AOFLAGGER_LIB)
  unset(AOFLAGGER_INCLUDE_DIR)
  if(AOFlagger_FIND_REQUIRED)
    message(FATAL_ERROR "AOFlagger required but not found.")
  else()
    message(INFO "AOFlagger not found, but not required.")
  endif()
else()
  message(STATUS "  AOFlagger include dir: ${AOFLAGGER_INCLUDE_DIR}")
  message(STATUS "  AOFlagger lib: ${AOFLAGGER_LIB}")
endif()
