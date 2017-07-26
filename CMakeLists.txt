# $Id$

lofar_add_package(DPPP)
lofar_add_package(TestDynDPPP)
lofar_add_package(PythonDPPP)
lofar_add_package(DPPP_AOFlag)
lofar_add_package(SPW_Combine SPWCombine)

lofar_find_package(Armadillo)
if(${ARMADILLO_FOUND})
  lofar_add_package(DPPP_DDECal)
else()
  message(WARNING "Armadillo was not found, NOT building DPPP_DDECal")
endif()
