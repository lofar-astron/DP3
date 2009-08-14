lib_LTLIBRARIES         = libdppp.la

libdppp_la_SOURCES = Package__Version.cc  \
		FlaggerStatistics.cc \
		MsInfo.cc \
		MsFile.cc \
		ComplexMedianFlagger.cc \
		FrequencyFlagger.cc \
		ComplexMedianFlagger2.cc \
		MADFlagger.cc \
		BandpassCorrector.cc \
		DataSquasher.cc \
		DataBuffer.cc \
		TimeBuffer.cc \
		RunDetails.cc \
		ProgressMeter.cc \
		Pipeline.cc


bin_PROGRAMS = versiondppp IDPPP

IDPPP_SOURCES	= IDPPP.cc PipelineProcessControl.cc
IDPPP_LDADD		= libdppp.la
IDPPP_DEPENDENCIES	= libdppp.la $(LOFAR_DEPEND)

pythondir = $(bindir)
dist_python_SCRIPTS =


versiondppp_SOURCES      = versiondppp.cc
versiondppp_LDADD        = libdppp.la
versiondppp_DEPENDENCIES = libdppp.la $(LOFAR_DEPEND)

include $(top_srcdir)/Makefile.common
