# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#[=======================================================================[.rst:
DetermineTargetCPU
------------------

Determine the target CPU used for the current build.

By default, the build system will generate code that is optimized for the
current architecture, by using the compile option ``-march=native``. This
will let the compiler determine the target CPU. The user many also manually
select a specific target CPU by setting the variable ``TARGET_CPU``. In order
to force the compiler to generate portable code, the user can set the option
``PORTABLE`` to ``TRUE``. Note that ``PORTABLE`` and ``TARGET_CPU`` are
mutually exclusive.

This modules sets the following ``INTERNAL`` variables:

::

  KNOWN_TARGET_CPUS - List of target CPUs known by both GCC and Clang
  IDENTIFIED_TARGET_CPU - Target CPU as identified by the compiler

#]=======================================================================]

# List of target CPUs known by both GCC and Clang This list was produced as
# follows (note: requires llc to be installed): comm -12 \ <(g++ -march=foo -E -
# < /dev/null |& grep '^cc1: note' | \ sed -nE 's,^.*: *([^;]*).*$,\1,p' | tr '
# ' '\n' | sort -u) \ <(llc -mattr=help |& grep processor. | awk '{print $1}' |
# sort -u)
set(KNOWN_TARGET_CPUS
    amdfam10
    athlon64
    athlon64-sse3
    athlon-fx
    atom
    barcelona
    bdver1
    bdver2
    bdver3
    bdver4
    bonnell
    broadwell
    btver1
    btver2
    core2
    core-avx2
    core-avx-i
    corei7
    corei7-avx
    haswell
    ivybridge
    k8
    k8-sse3
    knl
    nehalem
    nocona
    opteron
    opteron-sse3
    sandybridge
    silvermont
    skylake
    skylake-avx512
    slm
    westmere
    x86-64
    znver1
    CACHE INTERNAL "Known target CPUs")

if(NOT PORTABLE)
  get_directory_property(_compile_options COMPILE_OPTIONS)
  string(REPLACE ";" " " _compile_options "${_compile_options}")
  execute_process(
    COMMAND
      bash -c
      # Executed command is printed on stderr; we can discard stdout
      "echo | ${CMAKE_CXX_COMPILER} ${_compile_options} -E -v - >/dev/null"
    ERROR_VARIABLE _command
    RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(WARNING "${CMAKE_CXX_COMPILER_ID} compiler failed to identify "
                    "target CPU '${TARGET_CPU}'")
  else()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      execute_process(
        COMMAND
          bash -c
          "echo '${_command}' | sed -nE '/cc1/s/^.*-march=([^ ]+).*$/\\1/p'"
        OUTPUT_VARIABLE _target_cpu
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      execute_process(
        COMMAND
          bash -c
          "echo '${_command}' | sed -nE '/cc1/s/^.*-target-cpu ([^ ]+).*$/\\1/p'"
        OUTPUT_VARIABLE _target_cpu
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
      message(WARNING "Don't know how to let ${CMAKE_CXX_COMPILER_ID} "
                      "compiler identify target CPU")
    endif()
  endif()
endif()

if(DEFINED _target_cpu)
  set(IDENTIFIED_TARGET_CPU
      ${_target_cpu}
      CACHE INTERNAL "")
else()
  unset(IDENTIFIED_TARGET_CPU CACHE)
endif()

if(DEFINED IDENTIFIED_TARGET_CPU AND NOT IDENTIFIED_TARGET_CPU IN_LIST
                                     KNOWN_TARGET_CPUS)
  message(
    AUTHOR_WARNING
      "'${IDENTIFIED_TARGET_CPU}' is not in the list KNOWN_TARGET_CPUS. "
      "Please check if this list is still up-to-date")
endif()
