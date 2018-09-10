# Required for the ObjectFactory
include_directories(../../Common)

add_library(dppp_ddecal
  DDECal.cc Register.cc Stopwatch.cc
  KLFitter.cc DDECal.cc MultiDirSolver.cc Constraint.cc PiercePoint.cc
  ScreenConstraint.cc SmoothnessConstraint.cc
  TECConstraint.cc RotationConstraint.cc
  RotationAndDiagonalConstraint.cc
)
