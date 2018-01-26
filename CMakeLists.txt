# $Id$

lofar_add_package(DPPP)
lofar_add_package(TestDynDPPP)
lofar_add_package(PythonDPPP)
lofar_add_package(DPPP_AOFlag)
lofar_add_package(SPW_Combine SPWCombine)

lofar_find_package(Armadillo)
if(${ARMADILLO_FOUND})
  if(CMAKE_CXX_FLAGS MATCHES ".*\\+\\+11.*")
    lofar_add_package(DPPP_DDECal)
  else()
    message(WARNING "DPPP_DDECal will not build if you have no C++11 support, please build in directory gnucxx11_opt")
  endif()
else()
  message(WARNING "Armadillo was not found, NOT building DPPP_DDECal")
endif()
