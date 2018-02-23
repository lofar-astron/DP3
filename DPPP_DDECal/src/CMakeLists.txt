# $Id: CMakeLists.txt 30439 2014-11-19 15:04:34Z dijkema $

include(LofarPackageVersion)

lofar_add_library(dppp_ddecal
  Package__Version.cc
  DDECal.cc Register.cc Stopwatch.cc
  KLFitter.cc DDECal.cc MultiDirSolver.cc Constraint.cc PiercePoint.cc
  ScreenConstraint.cc TECConstraint.cc RotationConstraint.cc
  RotationAndDiagonalConstraint.cc
)

lofar_add_bin_program(versiondppp_ddecal versiondppp_ddecal.cc)
