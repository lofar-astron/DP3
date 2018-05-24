# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

lofar_package(DPPP_Interpolate 1.0 DEPENDS DPPP)

include(LofarFindPackage)
lofar_find_package(Casacore COMPONENTS casa ms tables REQUIRED)
lofar_find_package(Boost REQUIRED COMPONENTS date_time thread filesystem system python)

add_subdirectory(include/DPPP_Interpolate)
add_subdirectory(src)
add_subdirectory(test)
