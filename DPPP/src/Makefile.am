lib_LTLIBRARIES         = libdppp.la

libdppp_la_SOURCES = Package__Version.cc  \
		DPStep.cc DPInput.cc DPBuffer.cc AverageInfo.cc \
		ProgressMeter.cc FlagCounter.cc \
		UVWCalculator.cc \
		MSReader.cc MSWriter.cc MSUpdater.cc \
		Averager.cc MedFlagger.cc \
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
		Pipeline.cc


bin_PROGRAMS = versiondppp NDPPP IDPPP

versiondppp_SOURCES      = versiondppp.cc
versiondppp_LDADD        = libdppp.la
versiondppp_DEPENDENCIES = libdppp.la $(LOFAR_DEPEND)

NDPPP_SOURCES		= DPPP.cc
NDPPP_LDADD		= libdppp.la
NDPPP_DEPENDENCIES	= libdppp.la $(LOFAR_DEPEND)

IDPPP_SOURCES		= IDPPP.cc PipelineProcessControl.cc
IDPPP_LDADD		= libdppp.la
IDPPP_DEPENDENCIES	= libdppp.la $(LOFAR_DEPEND)

pythondir = $(bindir)
dist_python_SCRIPTS =

scriptdir = $(bindir)
dist_script_SCRIPTS = mwflagger mwflagger-part taqlflagger


include $(top_srcdir)/Makefile.common
