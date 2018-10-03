# Required for the ObjectFactory
include_directories(../Common)

include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(dppp_ddecal SHARED
  DDECal.cc Register.cc Stopwatch.cc
  KLFitter.cc DDECal.cc MultiDirSolver.cc Constraint.cc PiercePoint.cc
  ScreenConstraint.cc SmoothnessConstraint.cc
  TECConstraint.cc RotationConstraint.cc
  RotationAndDiagonalConstraint.cc
)
target_link_libraries(dppp_ddecal ${ARMADILLO_LIBRARIES})
set_target_properties(dppp_ddecal PROPERTIES VERSION 2)

install (TARGETS dppp_ddecal DESTINATION lib) 

#add_subdirectory(test)
