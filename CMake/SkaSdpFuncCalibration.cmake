# Provides the `sdp_func_calibration` target, which wraps the three-phase
# calibration aggregation of ska-sdp-func

include_guard(GLOBAL)

include(FetchContent)

# The aggregation is developed on a branch and is not part of any
# ska-sdp-func release yet, so this is pinned to a commit rather than a
# released tag. To move this on, check out the new upstream commit and
# reconfigure with -DSKA_SDP_FUNC_GIT_TAG=<sha>. Point this at a
# tag once the aggregation is part of a ska-sdp-func release.
#
# The default is a commit on master that moves the aggregation from
# OpenBLAS/LAPACK to Eigen (so it no longer needs a BLAS/LAPACKE backend).
set(SKA_SDP_FUNC_GIT_TAG
    "78f421b6f2d5b441e01fc9d347727737e5f73ec6"
    CACHE STRING "Git tag/commit providing the calibration aggregation")
set(SKA_SDP_FUNC_GIT_REPOSITORY
    "https://gitlab.com/ska-telescope/sdp/ska-sdp-func.git"
    CACHE STRING "Git repository providing the calibration aggregation")

message(STATUS "DP3: fetching ska-sdp-func for calibration aggregation")

# DP3 does not use any CUDA functionality, and a CUDA toolkit that happens to
# be installed can make the upstream build fail on its device code. EveryBeam
# does the same.
set(FIND_CUDA OFF CACHE INTERNAL "")

FetchContent_Declare(
  ska_sdp_func
  GIT_REPOSITORY ${SKA_SDP_FUNC_GIT_REPOSITORY}
  GIT_TAG ${SKA_SDP_FUNC_GIT_TAG}
  GIT_SHALLOW FALSE
  GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(ska_sdp_func)

# Upstream's compiler_options.cmake links the shared library with
# -static-libgcc -static-libstdc++ on GCC, which embeds a copy of libstdc++
# (including its __cxa_throw and transactional-memory clones) in
# libska_sdp_func. That copy then interposes the libstdc++.so.6 that the
# rest of DP3 uses, so an exception thrown in DP3 (for example an
# H5::Exception while opening an H5Parm) is dispatched through the wrong
# __cxa_throw and aborts. Hide the embedded libstdc++ symbols from the
# shared library's dynamic symbol table so they cannot interpose.
# --exclude-libs,ALL is a GNU-ld option; Apple's linker does not support it.
# AppleClang (CMAKE_CXX_COMPILER_ID == "AppleClang") is therefore excluded.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|^Clang$")
    target_link_options(ska_sdp_func PRIVATE -Wl,--exclude-libs,ALL)
    message(STATUS
            "DP3: hiding embedded libstdc++ symbols in the calibration "
            "aggregation")
endif()

# Upstream forces -stdlib=libc++ on its own sources when the compiler is
# Clang, while the rest of DP3 uses whatever the compiler uses by default.
# The two then have different C++ ABIs, and linking the aggregation library
# into DP3 fails with undefined references to std::__1:: symbols and a
# warning that libc++.so.1 was not found. Give the aggregation the same
# standard library as the rest of DP3: the target options come after the
# directory-wide flags, and the compiler uses the last -stdlib it is given.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_dp3_stdlib "libstdc++")
    if(CMAKE_CXX_FLAGS MATCHES "-stdlib=libc\\+\\+")
        set(_dp3_stdlib "libc++")
    endif()
    target_compile_options(ska_sdp_func PRIVATE "-stdlib=${_dp3_stdlib}")
    target_link_options(ska_sdp_func PRIVATE "-stdlib=${_dp3_stdlib}")
    message(STATUS
            "DP3: building the calibration aggregation with ${_dp3_stdlib}, "
            "to match the rest of DP3")
endif()

# Limit the number of symbols to link so that the system library do
# not clashs
get_target_property(_sdp_func_type ska_sdp_func TYPE)
if(_sdp_func_type STREQUAL "SHARED_LIBRARY"
   AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|^Clang$")
  target_link_options(ska_sdp_func INTERFACE -Wl,--exclude-libs,ALL)
endif()
