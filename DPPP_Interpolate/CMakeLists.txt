add_library(dppp_interpolate SHARED
  Interpolate.cc
)

set_target_properties(dppp_interpolate PROPERTIES VERSION 2)

install (TARGETS dppp_interpolate DESTINATION lib) 

#add_subdirectory(test)

