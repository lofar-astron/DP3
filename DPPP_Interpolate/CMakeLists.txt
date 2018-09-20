add_library(dppp_interpolate SHARED
  Interpolate.cc
)

install (TARGETS dppp_interpolate DESTINATION lib) 

#add_subdirectory(test)

