#[=======================================================================[.rst:
DetermineEnabledAVXFeatures
---------------------------

Determine which AVX-features are enabled by the current compile options.

By default, the build system will generate code that is optimized for the
current architecture, by using the compile option `-march=native`. This
module determines which AVX-features are enabled by the compiler. For each
enabled feature, a boolean variable ``USE_<AVX-feature>`` is written into
the cache. The user can force the generation of portable code by setting
the option ``PORTABLE`` to ``TRUE``. In that case, all cached 
``USE_<AVX-feature>`` variables will be removed from the cache.
#]=======================================================================]

get_cmake_property(_cache_variables CACHE_VARIABLES)
foreach(_var ${_cache_variables})
  if(_var MATCHES "^USE_AVX")
    if(PORTABLE)
      # We're building portable code, remove cached USE_AVX* variable
      unset(${_var} CACHE)
    elseif(${_var})
      # AVX-features were already determined or set explicitly. Bail out
      return()
    endif()
  endif()
endforeach()
get_directory_property(_compile_options COMPILE_OPTIONS)
string(REPLACE ";" " " _compile_options "${_compile_options}")
execute_process(COMMAND bash -c 
  "echo | ${CMAKE_CXX_COMPILER} ${_compile_options} -E -v - 2>&1 | grep cc1"
  OUTPUT_VARIABLE _compile_flags)
# Search for enabled AVX-features: GCC uses -m<feature> Clang +<feature>.
# Strip off the "-m" or "+" prefix to get the bare feature.
set(_pattern "(-m|\\+)(avx[0-9a-z]*)")
string(REGEX MATCHALL "${_pattern}" _avx_features "${_compile_flags}")
string(REGEX REPLACE "${_pattern}" "\\2" _avx_features "${_avx_features}")
foreach(_feature ${_avx_features})
  string(TOUPPER ${_feature} _FEATURE)
  set(USE_${_FEATURE} ON CACHE BOOL "Use ${_FEATURE}" FORCE)
  mark_as_advanced(USE_${_FEATURE})
endforeach()

