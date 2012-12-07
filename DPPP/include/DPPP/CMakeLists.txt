# $Id$

# List of header files that will be installed.
set(inst_HEADERS
  DPRun.h DPStep.h DPInput.h DPBuffer.h DPInfo.h
  DPLogger.h ProgressMeter.h FlagCounter.h
  UVWCalculator.h BaselineSelection.h
  MSReader.h MSWriter.h MSUpdater.h Counter.h
  Averager.h MedFlagger.h PreFlagger.h UVWFlagger.h AORFlagger.h
  StationAdder.h Filter.h
  PhaseShift.h Demixer.h
  Cursor.h CursorUtilCasa.h Position.h Stokes.h SourceDBUtil.h
  Apply.h EstimateMixed.h Simulate.h Simulator.h SubtractMixed.h Baseline.h
  ModelComponent.h PointSource.h GaussianSource.h Patch.h
  ModelComponentVisitor.h
  FlaggerStatistics.h)

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/${PACKAGE_NAME})

# Install header files.
install(FILES ${inst_HEADERS} DESTINATION include/${PACKAGE_NAME})
