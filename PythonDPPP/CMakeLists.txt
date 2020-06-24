# If Boost-Python3 found, build this package. Otherwise give a warning.
if(${Boost_PYTHON${BOOSTPY}_FOUND})

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
  find_package(PythonLibs 3 REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})

  if(${PYTHON_VERSION_MAJOR} EQUAL 2)
    find_package(Casacore COMPONENTS python)
  else()
    find_package(Casacore COMPONENTS python3)
  endif()
  
  find_package(Boost REQUIRED COMPONENTS python${BOOSTPY})
  include_directories(${Boost_INCLUDE_DIR})
  # By default, Boost-Python has auto_ptr support, which triggers a compiler
  # warning, since auto_ptr is deprecated -> Disable auto_ptr support.
  add_definitions(-DBOOST_NO_AUTO_PTR)

#  add_subdirectory(test)
else()
  message (WARNING, " Boost-Python not found; PythonDPPP will not be built.")
endif()
