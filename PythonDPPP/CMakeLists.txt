
pybind11_add_module(pydppp MODULE pydppp.cc)

add_library(PythonStep INTERFACE)

target_sources(PythonStep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/PyDPStep.cc)

target_link_libraries(PythonStep INTERFACE pybind11::embed)

add_executable(testpydpstep testmydpstep.cc
  ${DPPP_OBJECT}
  ${PARMDB_OBJECT}
  ${BLOB_OBJECT}
  ${COMMON_OBJECT}
  ${AOFLAGGERSTEP_OBJECT}
  ${DDECAL_OBJECT}
  ${IDGPREDICT_OBJECT}
)
target_link_libraries(testpydpstep PRIVATE PythonStep)
set_target_properties(testpydpstep PROPERTIES ENABLE_EXPORTS ON)


target_link_libraries(testpydpstep PRIVATE
  ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}
  ${CASACORE_LIBRARIES}
  ${LOFAR_STATION_RESPONSE_LIB}
  ${CFITSIO_LIBRARY}
  ${EXTRA_LIBRARIES})
