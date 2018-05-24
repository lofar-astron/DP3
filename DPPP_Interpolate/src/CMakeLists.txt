include(LofarPackageVersion)

lofar_add_library(dppp_interpolate
  Interpolate.cc ../../../../LCS/Common/src/LofarLogCout.cc
)
