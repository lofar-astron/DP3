# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

lofar_package(DPPP_AOFlag 1.0 DEPENDS DPPP)

include(LofarFindPackage)
lofar_find_package(AOFlagger REQUIRED)
lofar_find_package(Casacore COMPONENTS casa ms tables REQUIRED)
lofar_find_package(Boost REQUIRED COMPONENTS date_time thread filesystem system)
#lofar_find_package(GSL)
lofar_find_package(LibXml2 REQUIRED)
lofar_find_package(PNG REQUIRED)
lofar_find_package(FFTW3 REQUIRED)
lofar_find_package(CFITSIO REQUIRED)
#if(CMAKE_SYSTEM_NAME MATCHES "Linux")
#  find_library(RT_LIBRARY rt)
#  list(APPEND LOFAR_EXTRA_LIBRARIES ${RT_LIBRARY})
#endif(CMAKE_SYSTEM_NAME MATCHES "Linux")

add_subdirectory(include/DPPP_AOFlag)
add_subdirectory(src)
add_subdirectory(test)
