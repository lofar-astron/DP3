# Run script to find current directory (needed for getting a tarball from
# the source directory).
configure_file(findenv.run_tmpl findenv.run_script)

include(LofarCTest)

lofar_add_test(tDDECal)
