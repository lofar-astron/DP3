# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

lofar_package(DPPP_DDECal 1.0 DEPENDS DPPP)

include(LofarFindPackage)
lofar_find_package(Casacore COMPONENTS casa ms tables REQUIRED)
lofar_find_package(Armadillo REQUIRED)

add_subdirectory(include/DPPP_DDECal)
add_subdirectory(src)
add_subdirectory(test)
