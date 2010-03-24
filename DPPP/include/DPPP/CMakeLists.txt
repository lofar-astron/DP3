# $Id$

# List of header files that will be installed.
set(inst_HEADERS
  DPRun.h DPStep.h DPInput.h DPBuffer.h AverageInfo.h
  ProgressMeter.h FlagCounter.h
  UVWCalculator.h
  MSReader.h MSWriter.h MSUpdater.h Counter.h
  Averager.h MedFlagger.h PreFlagger.h UVWFlagger.h
  FlaggerStatistics.h
  MsInfo.h
  MsFile.h
  ComplexMedianFlagger.h
  FrequencyFlagger.h
  ComplexMedianFlagger2.h
  MADFlagger.h
  BandpassCorrector.h
  DataSquasher.h
  DataBuffer.h
  TimeBuffer.h
  RunDetails.h
  Pipeline.h)

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/${PACKAGE_NAME})

# Install header files.
install(FILES ${inst_HEADERS} DESTINATION include/${PACKAGE_NAME})
