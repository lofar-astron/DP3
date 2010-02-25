check_PROGRAMS        = tMedian tMirror tAverager tMedFlagger tNDPPP \
			test_flaggers

# programs to run through supplied checktools
CHECKTOOLPROGS        = $(CHECK_PROGRAMS)
#ENDCHECKTOOLPROGS

# scripts used to run tests
TESTSCRIPTS           = tMedian.sh tMirror.sh tAverager.sh tMedFlagger.sh \
			tNDPPP.sh \
			test_flaggers.sh tmwflagger.sh

# scripts and possible programs to run tests
TESTS                 = $(TESTSCRIPTS)

XFAIL_TESTS           =

# all files (.run, .stdout, .in, .log_prop, etc.) needed to run tests
EXTRA_DIST            = $(TESTSCRIPTS) \
			tmwflagger.run \
			tmwflagger.in_cd tmwflagger.in_vd \
			tmwflagger.parset tmwflagger.stdout

# Lines to build a test program testprg
tMedian_SOURCES                = tMedian.cc
tMedian_LDADD                  = ../src/libdppp.la
tMedian_DEPENDENCIES           = ../src/libdppp.la $(LOFAR_DEPEND)

tMirror_SOURCES                = tMirror.cc
tMirror_LDADD                  = ../src/libdppp.la
tMirror_DEPENDENCIES           = ../src/libdppp.la $(LOFAR_DEPEND)

tAverager_SOURCES              = tAverager.cc
tAverager_LDADD                = ../src/libdppp.la
tAverager_DEPENDENCIES         = ../src/libdppp.la $(LOFAR_DEPEND)

tMedFlagger_SOURCES            = tMedFlagger.cc
tMedFlagger_LDADD              = ../src/libdppp.la
tMedFlagger_DEPENDENCIES       = ../src/libdppp.la $(LOFAR_DEPEND)

tNDPPP_SOURCES                 = tNDPPP.cc
tNDPPP_LDADD                   = ../src/libdppp.la
tNDPPP_DEPENDENCIES	       = ../src/libdppp.la $(LOFAR_DEPEND)

test_flaggers_SOURCES                = test_flaggers.cc
test_flaggers_LDADD                  = ../src/libdppp.la
test_flaggers_DEPENDENCIES           = ../src/libdppp.la $(LOFAR_DEPEND)


include $(top_srcdir)/Makefile.common
