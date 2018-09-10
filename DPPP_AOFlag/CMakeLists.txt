# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

find_package(Casacore COMPONENTS casa ms tables REQUIRED)
find_package(Boost REQUIRED COMPONENTS date_time thread filesystem system python)
# AOFlagger depends on Python 2.7, see aoflagger CMake
find_package(Python 2.7 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(PNG REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_package(CFITSIO REQUIRED)

find_library(AOFLAGGER_LIB aoflagger)
find_path(AOFLAGGER_INCLUDE_DIR NAMES aoflagger.h)
include_directories(${AOFLAGGER_INCLUDE_DIR})

# Copied from AOFlagger CMakeLists to link in GTKMM if AOFlagger was built with it
find_package(PkgConfig)
pkg_check_modules(GTKMM gtkmm-3.0>=3.0.0)
pkg_check_modules(SIGCXX sigc++-2.0)

if(GTKMM_FOUND)
  set(LOFAR_EXTRA_LIBRARIES ${LOFAR_EXTRA_LIBRARIES} ${GTKMM_LIBRARIES} ${GLIBMM_LIBRARIES})
endif(GTKMM_FOUND)
# End check for GTKMM

add_subdirectory(src)
#add_subdirectory(test)
