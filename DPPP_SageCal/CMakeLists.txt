# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

lofar_package(DPPP_SageCal 1.0 DEPENDS DPPP)

include(LofarFindPackage)
lofar_find_package(SageCal REQUIRED)
lofar_find_package(BLAS REQUIRED)
lofar_find_package(Casacore COMPONENTS casa ms tables REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB_PKG glib-2.0)

include_directories(${SageCal_INCLUDE_DIR})
include_directories(${SageCal_INCLUDE_DIR}/sagecal)
include_directories(${GLIB_PKG_INCLUDE_DIRS})

set(${PACKAGE_NAME}_LINK_LIBRARIES ${${PACKAGE_NAME}_LINK_LIBRARIES} ${SageCal_LIBRARY} ${Dirac_LIBRARY} ${SageCalData_LIBRARY} ${GLIB_PKG_LIBRARIES})

add_subdirectory(include/DPPP_SageCal)
add_subdirectory(src)
