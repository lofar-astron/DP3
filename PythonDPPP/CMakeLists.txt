# Try to find Boost-Python.
# If found, build this package. Otherwise give a warning.
FIND_PATH(BOOST_PYTHON_FOUND "boost/python.hpp")

if(BOOST_PYTHON_FOUND)

  include(PythonInstall)

  add_library(_pythondppp MODULE PythonDPPP.cc)
  set_target_properties(_pythondppp PROPERTIES 
    SOVERSION ${DPPP_GLOBAL_SOVERSION}
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_BUILD_DIR}/lofar/pythondppp)

# This is a quick-and-dirty fix to install the Python binding module in the
# right place. It will now be installed twice, because lofar_add_library()
# will install it in $prefix/$libdir
  install(TARGETS _pythondppp
    DESTINATION ${PYTHON_INSTALL_DIR}/lofar/pythondppp)

# We only need this library because dependent packages need to link against
# this package's version info.
#add_library(dppp_pythondppp
#  PythonStep.cc
#  DPStepBase.cc
#)

  # Install Python modules
  python_install(__init__.py DESTINATION lofar/pythondppp)

	find_package(PythonInterp REQUIRED)
  find_package(PythonLibs 2.6 REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})

  find_package(Casacore COMPONENTS python)
  
  find_package(Boost REQUIRED COMPONENTS python)
  include_directories(${Boost_INCLUDE_DIR})

#  add_subdirectory(test)
else()
  message (WARNING, " Boost-Python not found; PythonDPPP will not be built.")
endif()
