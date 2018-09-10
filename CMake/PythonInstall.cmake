# - Install Python source files.
#  python_install(source1..sourceN DESTINATION install_dir)
# Install Python source files and byte-compile them in the directory
# ${PYTHON_INSTALL_DIR}/${install_dir}.

# Copyright (C) 2008-2009
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, softwaresupport@astron.nl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# $Id: PythonInstall.cmake 32905 2015-11-17 15:31:54Z schaap $

# Search for the Python interpreter.
find_package(PythonInterp)

# Derive the Python site-packages installation directory and build directory.
if(PYTHON_EXECUTABLE)
  set(_cmd
    "from distutils.sysconfig import get_python_lib"
    "from os.path import join"
    "print(join(
       get_python_lib(plat_specific=True, standard_lib=True, prefix=''), 
       'site-packages'))"
  )
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" "-c" "${_cmd}"
    OUTPUT_VARIABLE _pydir
    ERROR_VARIABLE _pyerr
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_pyerr)
    message(FATAL_ERROR "Python command failed:\n${_pyerr}")
  endif(_pyerr)
  
  if(NOT DEFINED PYTHON_BUILD_DIR)
    set(_PRINT_PYTHON_DIRS TRUE)
  endif()
  
  set(PYTHON_BUILD_DIR "${CMAKE_BINARY_DIR}/${_pydir}" CACHE PATH 
    "Build directory for Python extensions" FORCE)
  set(PYTHON_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${_pydir}" CACHE PATH 
    "Installation directory for Python extensions" FORCE)

  if(_PRINT_PYTHON_DIRS)
    message(STATUS "Build directory for Python extensions:        ${PYTHON_BUILD_DIR}")
    message(STATUS "Installation directory for Python extensions: ${PYTHON_INSTALL_DIR}")
  endif()
endif(PYTHON_EXECUTABLE)


#
# macro python_install
#
macro(python_install)

  # Precondition check.
  if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "python_install: Python interpreter not available")
  endif(NOT PYTHON_EXECUTABLE)

  # Parse arguments.
  # apart from the python files list, there are two additional arguments
  # DESTINATION (required), where to put the py files (relative to python lib dir)
  # EXECUTABLE (optional), makes the py files executable
  string(REGEX REPLACE ";?DESTINATION.*" "" _py_files "${ARGN}")
  string(REGEX REPLACE ";?EXECUTABLE.*" "" _py_files "${_py_files}")
  string(REGEX MATCH "DESTINATION;.*" _dest_dir "${ARGN}")
  string(REGEX REPLACE "^DESTINATION;" "" _dest_dir "${_dest_dir}")
  string(REGEX REPLACE ";?EXECUTABLE.*" "" _dest_dir "${_dest_dir}")
  string(REGEX MATCH "EXECUTABLE;" _executable "${ARGN}")

  #check if optional argument EXECUTABLE is set
  #if so, then install the _py_files as EXECUTABLE type (executable)
  #else as normal files (not executable)
  if("${_executable}" STRGREATER "")
    set(INSTALL_TYPE PROGRAMS)
  else()
    set(INSTALL_TYPE FILES)
  endif("${_executable}" STRGREATER "")

  if(_py_files MATCHES "^$")
    message(FATAL_ERROR "python_install: no sources files specified")
  endif(_py_files MATCHES "^$")
  if(_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";")
    message(FATAL_ERROR "python_install: destination directory invalid")
  endif(_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";")

  # Set python package build/install directory.
  set(_inst_dir "${PYTHON_INSTALL_DIR}/${_dest_dir}")
  set(_build_dir "${PYTHON_BUILD_DIR}/${_dest_dir}")

  # Install and byte-compile each Python file.
  foreach(_py ${_py_files})
    get_filename_component(_py_path ${_py} PATH)
    get_filename_component(_py_abs ${_py} ABSOLUTE)
    
    # check if _py is a path in CMAKE_BINARY_DIR. If so, then it is most likely a configured_file. 
    # then strip the CMAKE_CURRENT_BINARY_DIR prefix.
    if(${_py} MATCHES "^(${CMAKE_CURRENT_BINARY_DIR})")
      string(REGEX REPLACE "^(${CMAKE_CURRENT_BINARY_DIR}/)" "" _py "${_py}")
      get_filename_component(_py_path ${_py} PATH)
    endif()

    # Create a symlink to each Python file; needed to mimic install tree.
    file(MAKE_DIRECTORY ${_build_dir}/${_py_path})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${_py_abs} ${_build_dir}/${_py})
    install(${INSTALL_TYPE} ${_py_abs} DESTINATION ${_inst_dir}/${_py_path})
    if(USE_PYTHON_COMPILATION)
      set(_py_code
        "import py_compile, os"
        "destdir = os.environ.get('DESTDIR','')"
        "print('-- Byte-compiling: %s${_inst_dir}/${_py}' % destdir)"
        "py_compile.compile('%s${DESTDIR}${_inst_dir}/${_py}' % destdir, doraise=True)")
      install(CODE 
        "execute_process(COMMAND ${PYTHON_EXECUTABLE} -c \"${_py_code}\"
                       RESULT_VARIABLE _result)
       if(NOT _result EQUAL 0)
         message(FATAL_ERROR \"Byte-compilation FAILED: \$ENV{DESTDIR}${_inst_dir}/${_py}\")
       endif(NOT _result EQUAL 0)")
    endif(USE_PYTHON_COMPILATION)
  endforeach(_py ${_py_files})

  # Make sure that there's a __init__.py file in each build/install directory.
  string(REGEX REPLACE "/" ";" _dir_list ${_dest_dir})
  set(_init_dir)
  foreach(_dir ${_dir_list})
    set(_init_dir "${_init_dir}/${_dir}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E touch
      "${PYTHON_BUILD_DIR}${_init_dir}/__init__.py")
    install(CODE 
      "execute_process(COMMAND ${CMAKE_COMMAND} -E touch 
        \"\$ENV{DESTDIR}${PYTHON_INSTALL_DIR}${_init_dir}/__init__.py\")")
  endforeach(_dir ${_dir_list})

endmacro(python_install)
