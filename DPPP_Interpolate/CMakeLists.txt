add_library(dppp_interpolate SHARED
  Interpolate.cc
)

set_target_properties(dppp_interpolate PROPERTIES VERSION ${DPPP_GLOBAL_SOVERSION})

install (TARGETS dppp_interpolate DESTINATION lib) 

#add_subdirectory(test)

