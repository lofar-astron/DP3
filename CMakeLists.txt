cmake_minimum_required(VERSION 2.8.8)

project(DP3)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

add_subdirectory(DPPP)
#add_subdirectory(TestDynDPPP)
#add_subdirectory(PythonDPPP)
#add_subdirectory(DPPP_AOFlag)
#add_subdirectory(SPW_Combine SPWCombine)
#add_subdirectory(DPPP_Interpolate)

find_package(Armadillo)
if(${ARMADILLO_FOUND})
#  add_subdirectory(DPPP_DDECal)
else()
  message(WARNING "Armadillo was not found, NOT building DPPP_DDECal")
endif()
