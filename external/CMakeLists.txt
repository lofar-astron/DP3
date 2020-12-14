find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during configure" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        message(STATUS "Syncing submodules")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule sync --recursive
                        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "Syncing submodules did not succeed" 
            "command '${GIT_EXECUTABLE} submodule sync --recursive' failed with exit code ${GIT_SUBMOD_RESULT}")
        endif()
        message(STATUS "Updating submodules")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "Updating submodules did not succeed" 
            "command '${GIT_EXECUTABLE} submodule update --init --recursive' failed with exit code ${GIT_SUBMOD_RESULT}")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/pybind11/CMakeLists.txt")
    message(FATAL_ERROR "Submodules pybind11 was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules manually and try again.")
endif()

set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
add_subdirectory(pybind11)
