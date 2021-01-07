// DPRun.cc: Class to run steps like averaging and flagging on an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DPRun.h"

#include <boost/algorithm/string.hpp>

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "Averager.h"
#include "BDAAverager.h"
#include "ColumnReader.h"
#include "Counter.h"
#include "Demixer.h"
#include "DemixerNew.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Filter.h"
#include "GainCal.h"
#include "H5ParmPredict.h"
#include "Interpolate.h"
#include "MedFlagger.h"
#include "MSBDAWriter.h"
#include "MSReader.h"
#include "MSUpdater.h"
#include "MSWriter.h"
#include "MultiMSReader.h"
#include "PhaseShift.h"
#include "Predict.h"
#include "PreFlagger.h"
#include "ProgressMeter.h"
#include "ScaleData.h"
#include "SetBeam.h"
#include "Split.h"
#include "StationAdder.h"
#include "UVWFlagger.h"
#include "Upsample.h"

#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"

#include "../AOFlaggerStep/AOFlaggerStep.h"

#include "../DDECal/DDECal.h"

#include "../IDGPredict/IDGPredict.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/OS/DynLib.h>
#include <casacore/casa/Utilities/Regex.h>

namespace DP3 {
namespace DPPP {

// Initialize the statics.
std::map<std::string, DPRun::StepCtor*> DPRun::theirStepMap;

void DPRun::registerStepCtor(const std::string& type, StepCtor* func) {
  theirStepMap[type] = func;
}

DPRun::StepCtor* DPRun::findStepCtor(const std::string& type) {
  std::map<std::string, StepCtor*>::const_iterator iter =
      theirStepMap.find(type);
  if (iter != theirStepMap.end()) {
    return iter->second;
  }
  // Try to load the step from a dynamic library with that name
  // (in lowercase).
  // A dot can be used to have a specific library name (so multiple
  // steps can use the same shared library).
  std::string libname(type);
  boost::algorithm::to_lower(libname);
  string::size_type pos = libname.find_first_of(".");
  if (pos != string::npos) {
    libname = libname.substr(0, pos);
  }
  // Try to load and initialize the dynamic library.
  casacore::DynLib dl(libname, string("libdppp_"), "register_" + libname,
                      false);
  if (dl.getHandle()) {
    // See if registered now.
    iter = theirStepMap.find(type);
    if (iter != theirStepMap.end()) {
      return iter->second;
    }
  }
  throw Exception("Step type " + type +
                  " is unknown and no shared library lib" + libname +
                  " or libdppp_" + libname + " found in (DY)LD_LIBRARY_PATH");
}

void DPRun::execute(const string& parsetName, int argc, char* argv[]) {
  casacore::Timer timer;
  NSTimer nstimer;
  nstimer.start();
  ParameterSet parset;
  if (!parsetName.empty()) {
    parset.adoptFile(parsetName);
  }
  // Adopt possible parameters given at the command line.
  parset.adoptArgv(argc, argv);  ///< works fine if argc==0 and argv==0
  DPLogger::useLogger = parset.getBool("uselogger", false);
  bool showProgress = parset.getBool("showprogress", true);
  bool showTimings = parset.getBool("showtimings", true);
  // checkparset is an integer parameter now, but accept a bool as well
  // for backward compatibility.
  int checkparset = 0;
  try {
    checkparset = parset.getInt("checkparset", 0);
  } catch (...) {
    DPLOG_WARN_STR("Parameter checkparset should be an integer value");
    checkparset = parset.getBool("checkparset") ? 1 : 0;
  }

  bool showcounts = parset.getBool("showcounts", true);

  unsigned int numThreads = parset.getInt("numthreads", 0);

  // Create the steps, link them together
  DPInput::ShPtr firstStep = makeMainSteps(parset);

  DPStep::ShPtr step = firstStep;
  DPStep::ShPtr lastStep;
  while (step) {
    step = step->getNextStep();
  }

  // Call updateInfo()
  DPInfo dpInfo;
  if (numThreads > 0) {
    dpInfo.setNThreads(numThreads);
  }
  dpInfo = firstStep->setInfo(dpInfo);
  // Tell the reader if visibility data needs to be read.
  firstStep->setReadVisData(dpInfo.needVisData());

  // Show the steps.
  step = firstStep;
  while (step) {
    std::ostringstream os;
    step->show(os);
    DPLOG_INFO(os.str(), true);
    lastStep = step;
    step = step->getNextStep();
  }
  if (checkparset >= 0) {
    // Show unused parameters (might be misspelled).
    std::vector<std::string> unused = parset.unusedKeys();
    if (!unused.empty()) {
      DPLOG_WARN_STR(
          "\n*** WARNING: the following parset keywords were not used ***"
          << "\n             maybe they are misspelled"
          << "\n    " << unused << std::endl);
      if (checkparset != 0) throw Exception("Unused parset keywords found");
    }
  }
  // Process until the end.
  unsigned int ntodo = firstStep->getInfo().ntime();
  DPLOG_INFO_STR("Processing " << ntodo << " time slots ...");
  {
    ProgressMeter* progress = 0;
    if (showProgress) {
      progress = new ProgressMeter(0.0, ntodo, "NDPPP", "Time slots processed",
                                   "", "", true, 1);
    }
    double ndone = 0;
    if (showProgress && ntodo > 0) {
      progress->update(ndone, true);
    }
    DPBuffer buf;
    while (firstStep->process(buf)) {
      ++ndone;
      if (showProgress && ntodo > 0) {
        progress->update(ndone, true);
      }
    }
    delete progress;
  }
  // Finish the processing.
  DPLOG_INFO_STR("Finishing processing ...");
  firstStep->finish();
  // Give all steps the option to add something to the MS written.
  // It starts with the last step to get the name of the output MS,
  // but each step must first call its previous step before
  // it adds something itself.
  lastStep->addToMS("");

  // Show the counts where needed.
  if (showcounts) {
    step = firstStep;
    while (step) {
      std::ostringstream os;
      step->showCounts(os);
      DPLOG_INFO(os.str(), true);
      step = step->getNextStep();
    }
  }
  // Show the overall timer.
  nstimer.stop();
  double duration = nstimer.getElapsed();
  std::ostringstream ostr;
  ostr << std::endl;
  // Output special line for pipeline use.
  if (DPLogger::useLogger) {
    ostr << "Start timer output" << std::endl;
  }
  timer.show(ostr, "Total NDPPP time");
  DPLOG_INFO(ostr.str(), true);
  if (showTimings) {
    // Show the timings per step.
    step = firstStep;
    while (step) {
      std::ostringstream os;
      step->showTimings(os, duration);
      if (!os.str().empty()) {
        DPLOG_INFO(os.str(), true);
      }
      step = step->getNextStep();
    }
  }
  if (DPLogger::useLogger) {
    ostr << "End timer output\n";
  }
  // The destructors are called automatically at this point.
}

DPInput::ShPtr DPRun::makeMainSteps(const ParameterSet& parset) {
  DPInput::ShPtr inputStep = DPInput::CreateReader(parset, "");
  DPStep::ShPtr step = makeStepsFromParset(parset, "", inputStep.get(), false);
  if (step) inputStep->setNextStep(step);

  // Calculate needsOutputStep to be true if one of the steps changes
  // the visibility stream and therefore requires the output to be rewritten
  bool needsOutputStep = false;
  step = inputStep;
  while (step->getNextStep()) {
    step = step->getNextStep();
    needsOutputStep = needsOutputStep || step->modifiesData();
  }
  // step now points to the last step in the chain, which is required later on.

  // Add an output step if not explicitly added as the last step
  const bool endsWithOutputStep = dynamic_cast<MSWriter*>(step.get()) ||
                                  dynamic_cast<MSUpdater*>(step.get()) ||
                                  dynamic_cast<Split*>(step.get());

  if (!endsWithOutputStep) {
    const std::string msOutName = parset.getString("msout");
    if (needsOutputStep || !msOutName.empty()) {
      std::string msName = casacore::Path(inputStep->msName()).absoluteName();
      DPStep::ShPtr outputStep =
          makeOutputStep(inputStep.get(), parset, "msout.", msName,
                         step->outputs() == DPStep::MSType::BDA);
      step->setNextStep(outputStep);
      step = outputStep;
    }
  }

  // Add a null step, so the last step can use getNextStep->process().
  step->setNextStep(std::make_shared<NullStep>());
  return inputStep;
}

DPStep::ShPtr DPRun::makeStepsFromParset(const ParameterSet& parset,
                                         const std::string& prefix,
                                         DPInput* inputStep,
                                         bool terminateChain) {
  std::string msName;
  const std::vector<string> stepNames =
      parset.getStringVector(prefix + "steps");

  DPStep::ShPtr firstStep;
  DPStep::ShPtr lastStep;
  for (const std::string& stepName : stepNames) {
    std::string prefix(stepName + '.');

    // The alphabetic part of the name is the default step type.
    // This allows names like average1, out3.
    std::string defaultType = stepName;
    while (!defaultType.empty() && std::isdigit(defaultType.back())) {
      defaultType.resize(defaultType.size() - 1);
    }
    std::string type = parset.getString(prefix + "type", defaultType);
    boost::algorithm::to_lower(type);

    DPStep::MSType inputType = lastStep ? lastStep->outputs() : DPStep::REGULAR;
    DPStep::ShPtr step =
        makeSingleStep(type, inputStep, parset, prefix, msName, inputType);

    if (lastStep) {
      if (!step->accepts(lastStep->outputs())) {
        throw std::invalid_argument("Step " + type +
                                    " is incompatible with the input data.");
      }
      lastStep->setNextStep(step);
    }
    lastStep = step;

    if (!firstStep) {
      firstStep = step;
    }
  }

  if (terminateChain) {
    // Add a null step, so the last step can use getNextStep->process().
    lastStep->setNextStep(std::make_shared<NullStep>());
  }

  return firstStep;
}

DPStep::ShPtr DPRun::makeSingleStep(const std::string& type, DPInput* inputStep,
                                    const ParameterSet& parset,
                                    const std::string& prefix,
                                    std::string& msName,
                                    DPStep::MSType inputType) {
  if (type == "aoflagger" || type == "aoflag") {
    return std::make_shared<AOFlaggerStep>(inputStep, parset, prefix);
  } else if (type == "averager" || type == "average" || type == "squash") {
    return std::make_shared<Averager>(inputStep, parset, prefix);
  } else if (type == "bdaaverager") {
    return std::make_shared<BDAAverager>(*inputStep, parset, prefix);
  } else if (type == "madflagger" || type == "madflag") {
    return std::make_shared<MedFlagger>(inputStep, parset, prefix);
  } else if (type == "preflagger" || type == "preflag") {
    return std::make_shared<PreFlagger>(inputStep, parset, prefix);
  } else if (type == "uvwflagger" || type == "uvwflag") {
    return std::make_shared<UVWFlagger>(inputStep, parset, prefix);
  } else if (type == "columnreader") {
    return std::make_shared<ColumnReader>(*inputStep, parset, prefix);
  } else if (type == "counter" || type == "count") {
    return std::make_shared<Counter>(inputStep, parset, prefix);
  } else if (type == "phaseshifter" || type == "phaseshift") {
    return std::make_shared<PhaseShift>(inputStep, parset, prefix);
  } else if (type == "demixer" || type == "demix") {
    return std::make_shared<Demixer>(inputStep, parset, prefix);
  } else if (type == "smartdemixer" || type == "smartdemix") {
    return std::make_shared<DemixerNew>(inputStep, parset, prefix);
  } else if (type == "applybeam") {
    return std::make_shared<ApplyBeam>(inputStep, parset, prefix);
  } else if (type == "stationadder" || type == "stationadd") {
    return std::make_shared<StationAdder>(inputStep, parset, prefix);
  } else if (type == "scaledata") {
    return std::make_shared<ScaleData>(inputStep, parset, prefix);
  } else if (type == "setbeam") {
    return std::make_shared<SetBeam>(inputStep, parset, prefix);
  } else if (type == "filter") {
    return std::make_shared<Filter>(inputStep, parset, prefix);
  } else if (type == "applycal" || type == "correct") {
    return std::make_shared<ApplyCal>(inputStep, parset, prefix);
  } else if (type == "predict") {
    return std::make_shared<Predict>(inputStep, parset, prefix);
  } else if (type == "idgpredict") {
    return std::make_shared<IDGPredict>(*inputStep, parset, prefix);
  } else if (type == "h5parmpredict") {
    return std::make_shared<H5ParmPredict>(inputStep, parset, prefix);
  } else if (type == "gaincal" || type == "calibrate") {
    return std::make_shared<GainCal>(inputStep, parset, prefix);
  } else if (type == "upsample") {
    return std::make_shared<Upsample>(inputStep, parset, prefix);
  } else if (type == "split" || type == "explode") {
    return std::make_shared<Split>(inputStep, parset, prefix);
  } else if (type == "ddecal") {
    return std::make_shared<DDECal>(inputStep, parset, prefix);
  } else if (type == "interpolate") {
    return std::make_shared<Interpolate>(inputStep, parset, prefix);
  } else if (type == "out" || type == "output" || type == "msout") {
    if (msName.empty())
      msName = casacore::Path(inputStep->msName()).absoluteName();
    return makeOutputStep(inputStep, parset, prefix, msName,
                          inputType == DPStep::MSType::BDA);
  } else {
    // Maybe the step is defined in a dynamic library.
    return findStepCtor(type)(inputStep, parset, prefix);
  }
}

DPStep::ShPtr DPRun::makeOutputStep(DPInput* reader, const ParameterSet& parset,
                                    const std::string& prefix,
                                    std::string& currentMSName,
                                    const bool& isBDA) {
  DPStep::ShPtr step;
  std::string outName;
  bool doUpdate = false;
  if (prefix == "msout.") {
    // The last output step.
    outName = parset.getString("msout.name", "");
    if (outName.empty()) {
      outName = parset.getString("msout", "");
    }
  } else {
    // An intermediate output step.
    outName = parset.getString(prefix + "name");
  }

  // A name equal to . or the last name means an update of the last MS.
  if (outName.empty() || outName == ".") {
    outName = currentMSName;
    doUpdate = true;
  } else {
    casacore::Path pathOut(outName);
    if (currentMSName == std::string(pathOut.absoluteName())) {
      outName = currentMSName;
      doUpdate = true;
    }
  }
  if (isBDA) {
    if (doUpdate) {
      throw std::invalid_argument("No updater for BDA data.");
    } else {
      step = std::make_shared<MSBDAWriter>(reader, outName, parset, prefix);
      reader->setReadVisData(true);
    }
  } else {
    if (doUpdate) {
      // Create MSUpdater.
      // Take care the history is not written twice.
      // Note that if there is nothing to write, the updater won't do anything.
      step = std::make_shared<MSUpdater>(reader, outName, parset, prefix,
                                         outName != currentMSName);
    } else {
      step = std::make_shared<MSWriter>(reader, outName, parset, prefix);
      reader->setReadVisData(true);
    }
  }
  casacore::Path pathOut(outName);
  currentMSName = pathOut.absoluteName();
  return step;
}

}  // namespace DPPP
}  // namespace DP3
