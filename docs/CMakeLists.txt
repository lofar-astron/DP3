cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../CMake)

project(DP3_DOC NONE)

find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_LIST_DIR}/doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

find_package(Sphinx)
find_program(JSONSCHEMA2RST_EXECUTABLE
             NAMES jsonschema2rst
             DOC "Path to jsonschema2rst executable")

if(SPHINX_FOUND AND JSONSCHEMA2RST_EXECUTABLE)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_LIST_DIR})
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/)
  add_custom_command(OUTPUT inc_with_index COMMAND
                     ${JSONSCHEMA2RST_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/schemas ${CMAKE_CURRENT_SOURCE_DIR}/inc)
  add_custom_command(OUTPUT inc COMMAND
                     rm ${CMAKE_CURRENT_SOURCE_DIR}/inc/index.rst
                     DEPENDS inc_with_index)
  add_custom_target(userdoc
                    COMMAND
                    ${SPHINX_EXECUTABLE} -b html
                    ${SPHINX_SOURCE} ${SPHINX_BUILD}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    DEPENDS inc
                    COMMENT "Generating user documentation with Sphinx"
                    )
else()
  message("Sphinx (sphinx-build) and jsonschema2rst are necessary for building user docs")
endif(SPHINX_FOUND AND JSONSCHEMA2RST_EXECUTABLE)
