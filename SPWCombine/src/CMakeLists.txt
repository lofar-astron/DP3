# $Id$

include(LofarPackageVersion)

lofar_add_library(spwcombine
  Package__Version.cc)

lofar_add_bin_program(SPWCombine 
  SPWCombine.cc
  Combiner.cc
  CombinerProcessControl.cc)

lofar_add_bin_program(versionspw_combine versionspw_combine.cc)
