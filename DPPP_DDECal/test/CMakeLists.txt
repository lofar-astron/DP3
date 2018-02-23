# $Id: CMakeLists.txt 38047 2017-07-26 12:49:53Z dijkema $

# Run script to find current directory (needed for getting a tarball from
# the source directory).
configure_file(findenv.run_tmpl findenv.run_script)

include(LofarCTest)

lofar_add_test(tDDECal)
lofar_add_test(tRotationConstraint tRotationConstraint.cc)
