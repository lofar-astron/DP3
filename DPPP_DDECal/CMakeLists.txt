# $Id: CMakeLists.txt 27640 2013-12-04 08:02:49Z diepen $

find_package(Casacore COMPONENTS casa ms tables REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS date_time)

add_subdirectory(src)
#add_subdirectory(test)
