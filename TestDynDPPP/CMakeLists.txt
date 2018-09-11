add_library(dppp_testdyndppp SHARED
  TestDynStep.cc ${DPPP_OBJECT} ${PARMDB_OBJECT} ${BLOB_OBJECT} ${COMMON_OBJECT}
)

#add_subdirectory(test)
