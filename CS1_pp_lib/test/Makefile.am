check_PROGRAMS        = test_flaggers

# programs to run through supplied checktools
CHECKTOOLPROGS        = $(CHECK_PROGRAMS)
#ENDCHECKTOOLPROGS

# scripts used to run tests
TESTSCRIPTS           = test_flaggers.sh

# scripts and possible programs to run tests
TESTS                 = $(TESTSCRIPTS)

XFAIL_TESTS           =

# all files (.run, .stdout, .in, .log_prop, etc.) needed to run tests
EXTRA_DIST            = $(TESTSCRIPTS)

# Lines to build a test program testprg
test_flaggers_SOURCES                = test_flaggers.cc
test_flaggers_LDADD                  = ../src/libcs1_pp_lib.la
test_flaggers_DEPENDENCIES           = ../src/libcs1_pp_lib.la $(LOFAR_DEPEND)


include $(top_srcdir)/Makefile.common
