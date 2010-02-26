INSTHDRS = \
		DPRun.h DPStep.h DPInput.h DPBuffer.h AverageInfo.h \
		ProgressMeter.h FlagCounter.h \
		UVWCalculator.h \
		MSReader.h MSWriter.h MSUpdater.h \
		Averager.h MedFlagger.h PreFlagger.h \
		FlaggerStatistics.h \
		MsInfo.h \
		MsFile.h \
		ComplexMedianFlagger.h \
		FrequencyFlagger.h \
		ComplexMedianFlagger2.h \
		MADFlagger.h \
		BandpassCorrector.h \
		DataSquasher.h \
		DataBuffer.h \
		TimeBuffer.h \
		RunDetails.h \
		Pipeline.h 

NOINSTHDRS = PipelineProcessControl.h 

TCCHDRS =

pkginclude_HEADERS = Package__Version.h $(INSTHDRS) $(TCCHDRS)

noinst_HEADERS = $(NOINSTHDRS)

DOCHDRS = $(INSTHDRS) $(NOINSTHDRS)


include $(top_srcdir)/Makefile.common
