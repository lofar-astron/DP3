# FindNVML.cmake

if(${CUDA_VERSION_STRING} VERSION_LESS "9.1")
    string(CONCAT ERROR_MSG "--> ARCHER: Current CUDA version "
                         ${CUDA_VERSION_STRING}
                         " is too old. Must upgrade it to 9.1 or newer.")
    message(FATAL_ERROR ${ERROR_MSG})
endif()

# windows, including both 32-bit and 64-bit
if(WIN32)
    set(NVML_NAMES nvml)
    set(NVML_LIB_DIR "${CUDA_TOOLKIT_ROOT_DIR}/lib/x64")
    set(NVML_INCLUDE_DIR ${CUDA_INCLUDE_DIRS})

    # .lib import library full path
    find_file(NVML_LIB_PATH
              NO_DEFAULT_PATH
              NAMES nvml.lib
              PATHS ${NVML_LIB_DIR})

    # .dll full path
    find_file(NVML_DLL_PATH
              NO_DEFAULT_PATH
              NAMES nvml.dll
              PATHS "C:/Program Files/NVIDIA Corporation/NVSMI")
# linux
elseif(UNIX AND NOT APPLE)
    set(NVML_NAMES nvidia-ml)
    set(NVML_LIB_DIR "${CUDA_TOOLKIT_ROOT_DIR}/lib64/stubs" "${CUDA_TOOLKIT_ROOT_DIR}/lib/x86_64-linux-gnu")
    set(NVML_INCLUDE_DIR ${CUDA_INCLUDE_DIRS})

    find_library(NVML_LIB_PATH
                 NO_DEFAULT_PATH
                 NAMES ${NVML_NAMES}
                 PATHS ${NVML_LIB_DIR})
else()
    message(FATAL_ERROR "Unsupported platform.")
endif()

find_path(NVML_INCLUDE_PATH
          NO_DEFAULT_PATH
          NAMES nvml.h
          PATHS ${NVML_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NVML DEFAULT_MSG NVML_LIB_PATH NVML_INCLUDE_PATH)
