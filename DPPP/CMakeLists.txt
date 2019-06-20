include(PythonInstall)

add_executable(DPPP DPPP.cc
  ${DPPP_OBJECT}
  ${PARMDB_OBJECT}
  ${BLOB_OBJECT}
  ${COMMON_OBJECT}
  ${AOFLAGGERSTEP_OBJECT}
  ${DDECAL_OBJECT}
  ${IDGPREDICT_OBJECT})
target_link_libraries(DPPP
  ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}
  ${CASACORE_LIBRARIES}
  ${LOFAR_STATION_RESPONSE_LIB}
  ${CFITSIO_LIBRARY}
  ${EXTRA_LIBRARIES})

install(TARGETS DPPP DESTINATION bin)

# Data files
install(
  FILES share/LOFAR-LBA-default.rfis share/LOFAR-HBA-default.rfis
  DESTINATION share/rfistrategies)

# lofar_add_bin_scripts(taqlflagger)

# Python modules.
python_install(
  __init__.py
  DESTINATION lofar/dppp)
