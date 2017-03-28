# $Id$

lofar_package(DPPP 1.0 DEPENDS LofarStMan Common MS ParmDB StationResponse)

include(LofarFindPackage)
lofar_find_package(Casacore COMPONENTS casa ms tables REQUIRED)
lofar_find_package(HDF5 COMPONENTS CXX REQUIRED)

add_subdirectory(include/DPPP)
add_subdirectory(src)
add_subdirectory(share)
add_subdirectory(etc)
add_subdirectory(test)
