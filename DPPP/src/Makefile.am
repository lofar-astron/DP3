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
		RunDetails.cc \
		Pipeline.cc 


bin_PROGRAMS = versiondppp CS1_IDPPP

CS1_IDPPP_SOURCES	= IDPPP.cc PipelineProcessControl.cc
CS1_IDPPP_LDADD		= libdppp.la
CS1_IDPPP_DEPENDENCIES	= libdppp.la $(LOFAR_DEPEND)

pythondir = $(bindir)
dist_python_SCRIPTS = 


versiondppp_SOURCES      = versiondppp.cc
versiondppp_LDADD        = libdppp.la
versiondppp_DEPENDENCIES = libdppp.la $(LOFAR_DEPEND)

include $(top_srcdir)/Makefile.common
