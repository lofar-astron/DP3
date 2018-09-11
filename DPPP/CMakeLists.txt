include(PythonInstall)

# Required for the ObjectFactory
#include_directories(../Common)

add_executable(DPPP NDPPP.cc ${DPPP_OBJECT} ${PARMDB_OBJECT} ${BLOB_OBJECT} ${COMMON_OBJECT})
target_link_libraries(DPPP ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${CASACORE_LIBRARIES} ${LOFAR_STATION_RESPONSE_LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# lofar_add_bin_scripts(taqlflagger)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)

add_subdirectory(share)
#add_subdirectory(test)
