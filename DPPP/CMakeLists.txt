find_package(Casacore COMPONENTS casa ms tables REQUIRED)

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS filesystem system REQUIRED)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

add_subdirectory(src)
add_subdirectory(share)
#add_subdirectory(test)
