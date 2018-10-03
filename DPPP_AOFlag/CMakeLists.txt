find_package(Boost REQUIRED COMPONENTS date_time filesystem system python)
# AOFlagger depends on Python 2.7, see aoflagger CMake
find_package(PythonInterp REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
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
  set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} ${GTKMM_LIBRARIES} ${GLIBMM_LIBRARIES})
endif(GTKMM_FOUND)

add_library(dppp_aoflag SHARED
  AOFlaggerStep.cc
  Register.cc
)
target_link_libraries(dppp_aoflag
  ${CASACORE_LIBRARIES}
  ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_PYTHON_LIBRARY}
  ${EXTRA_LIBRARIES}
  ${AOFLAGGER_LIB})

install (TARGETS dppp_aoflag DESTINATION lib) 

#add_subdirectory(test)
