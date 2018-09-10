cmake_minimum_required(VERSION 2.8.8)

project(DP3)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/Station.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
	include_directories(${LOFAR_STATION_RESPONSE_DIR})
  set(LBEAM_LIBS ${LOFAR_STATION_RESPONSE_LIB})
  set(LBEAM_FILES lofar/lbeamevaluator.cpp lofar/lmspredicter.cpp aterms/lofarbeamterm.cpp)
  message("LOFAR beam library found.")
else()
  message(FATAL_ERROR "LOFAR beam library not found")
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

add_subdirectory(DPPP)
#add_subdirectory(TestDynDPPP)
#add_subdirectory(PythonDPPP)
add_subdirectory(DPPP_AOFlag)
#add_subdirectory(SPW_Combine SPWCombine)
#add_subdirectory(DPPP_Interpolate)

find_package(Armadillo)
if(${ARMADILLO_FOUND})
#  add_subdirectory(DPPP_DDECal)
else()
  message(WARNING "Armadillo was not found, NOT building DPPP_DDECal")
endif()
