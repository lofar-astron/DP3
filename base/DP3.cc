// DPRun.cc: Class to run steps like averaging and flagging on an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DP3.h"

#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "ProgressMeter.h"

#include <boost/algorithm/string.hpp>

#include "../steps/AOFlaggerStep.h"
#include "../steps/ApplyBeam.h"
#include "../steps/ApplyCal.h"
#include "../steps/Averager.h"
#include "../steps/BDAAverager.h"
#include "../steps/BDAExpander.h"
#include "../steps/BdaGroupPredict.h"
#include "../steps/ColumnReader.h"
#include "../steps/Counter.h"
#include "../steps/DDECal.h"
#include "../steps/BdaDdeCal.h"
#include "../steps/Demixer.h"
#include "../steps/DemixerNew.h"
#include "../steps/Filter.h"
#include "../steps/GainCal.h"
#include "../steps/H5ParmPredict.h"
#include "../steps/IDGPredict.h"
#include "../steps/Interpolate.h"
#include "../steps/MedFlagger.h"
#include "../steps/MSBDAWriter.h"
#include "../steps/MSReader.h"
#include "../steps/MSUpdater.h"
#include "../steps/MSWriter.h"
#include "../steps/MultiMSReader.h"
#include "../steps/PhaseShift.h"
#include "../steps/Predict.h"
#include "../steps/PreFlagger.h"
#include "../steps/ScaleData.h"
#include "../steps/SetBeam.h"
#include "../steps/Split.h"
#include "../steps/StationAdder.h"
#include "../steps/UVWFlagger.h"
#include "../steps/Upsample.h"

#include "../pythondp3/PyStep.h"

#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/OS/DynLib.h>
#include <casacore/casa/Utilities/Regex.h>

using dp3::steps::InputStep;
using dp3::steps::MSBDAWriter;
using dp3::steps::MSUpdater;
using dp3::steps::MSWriter;
using dp3::steps::NullStep;
using dp3::steps::Split;
using dp3::steps::Step;
using dp3::common::operator<<;

namespace dp3 {
namespace base {

// Initialize the statics.
std::map<std::string, DP3::StepCtor*> DP3::theirStepMap;

void DP3::registerStepCtor(const std::string& type, StepCtor* func) {
  theirStepMap[type] = func;
}

DP3::StepCtor* DP3::findStepCtor(const std::string& type) {
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

void DP3::execute(const string& parsetName, int argc, char* argv[]) {
  casacore::Timer timer;
  common::NSTimer nstimer;
  nstimer.start();
  common::ParameterSet parset;
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
  std::shared_ptr<InputStep> firstStep = makeMainSteps(parset);

  Step::ShPtr step = firstStep;
  Step::ShPtr lastStep;
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
  if (showProgress) {
    double ndone = 0;
    ProgressMeter progress(ndone, ntodo, "DP3", "Time slots processed", "", "",
                           true, 1);
    if (ntodo > 0) progress.update(ndone, true);
    DPBuffer buf;
    while (firstStep->process(buf)) {
      ++ndone;
      if (ntodo > 0) progress.update(ndone, true);
    }
  } else {
    DPBuffer buf;
    while (firstStep->process(buf)) {
    }
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

std::shared_ptr<InputStep> DP3::makeMainSteps(
    const common::ParameterSet& parset) {
  std::shared_ptr<InputStep> inputStep = InputStep::CreateReader(parset);
  Step::ShPtr step = makeStepsFromParset(parset, "", "steps", *inputStep, false,
                                         inputStep->outputs());
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
    const std::string msOutName = parset.getString(
        parset.isDefined("msout.name") ? "msout.name" : "msout");
    if (needsOutputStep || !msOutName.empty()) {
      std::string msName = casacore::Path(inputStep->msName()).absoluteName();
      Step::ShPtr outputStep = makeOutputStep(inputStep.get(), parset, "msout.",
                                              msName, step->outputs());
      step->setNextStep(outputStep);
      step = outputStep;
    }
  }

  // Add a null step, so the last step can use getNextStep->process().
  step->setNextStep(std::make_shared<NullStep>());
  return inputStep;
}

Step::ShPtr DP3::makeStepsFromParset(const common::ParameterSet& parset,
                                     const std::string& prefix,
                                     const std::string& step_names_key,
                                     InputStep& inputStep, bool terminateChain,
                                     Step::MsType initial_step_output) {
  std::string msName;
  const std::vector<string> stepNames =
      parset.getStringVector(prefix + step_names_key);

  Step::ShPtr firstStep;
  Step::ShPtr lastStep;
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

    Step::MsType inputType =
        lastStep ? lastStep->outputs() : initial_step_output;
    Step::ShPtr step =
        makeSingleStep(type, &inputStep, parset, prefix, msName, inputType);

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

std::shared_ptr<Step> DP3::makeSingleStep(const std::string& type,
                                          InputStep* inputStep,
                                          const common::ParameterSet& parset,
                                          const std::string& prefix,
                                          std::string& msName,
                                          Step::MsType inputType) {
  std::shared_ptr<Step> step;
  if (type == "aoflagger" || type == "aoflag") {
    step = std::make_shared<steps::AOFlaggerStep>(inputStep, parset, prefix);
  } else if (type == "averager" || type == "average" || type == "squash") {
    step = std::make_shared<steps::Averager>(*inputStep, parset, prefix);
  } else if (type == "bdaaverage" || type == "bdaaverager") {
    step = std::make_shared<steps::BDAAverager>(*inputStep, parset, prefix);
  } else if (type == "bdaexpander") {
    step = std::make_shared<steps::BDAExpander>(prefix);
  } else if (type == "madflagger" || type == "madflag") {
    step = std::make_shared<steps::MedFlagger>(inputStep, parset, prefix);
  } else if (type == "preflagger" || type == "preflag") {
    step = std::make_shared<steps::PreFlagger>(inputStep, parset, prefix);
  } else if (type == "uvwflagger" || type == "uvwflag") {
    step = std::make_shared<steps::UVWFlagger>(inputStep, parset, prefix);
  } else if (type == "columnreader") {
    step = std::make_shared<steps::ColumnReader>(*inputStep, parset, prefix);
  } else if (type == "counter" || type == "count") {
    step = std::make_shared<steps::Counter>(inputStep, parset, prefix);
  } else if (type == "phaseshifter" || type == "phaseshift") {
    step = std::make_shared<steps::PhaseShift>(inputStep, parset, prefix);
  } else if (type == "demixer" || type == "demix") {
    step = std::make_shared<steps::Demixer>(inputStep, parset, prefix);
  } else if (type == "smartdemixer" || type == "smartdemix") {
    step = std::make_shared<steps::DemixerNew>(inputStep, parset, prefix);
  } else if (type == "applybeam") {
    step = std::make_shared<steps::ApplyBeam>(inputStep, parset, prefix);
  } else if (type == "stationadder" || type == "stationadd") {
    step = std::make_shared<steps::StationAdder>(inputStep, parset, prefix);
  } else if (type == "scaledata") {
    step = std::make_shared<steps::ScaleData>(inputStep, parset, prefix,
                                              inputType);
  } else if (type == "setbeam") {
    step = std::make_shared<steps::SetBeam>(inputStep, parset, prefix);
  } else if (type == "filter") {
    step = std::make_shared<steps::Filter>(inputStep, parset, prefix);
  } else if (type == "applycal" || type == "correct") {
    step = std::make_shared<steps::ApplyCal>(inputStep, parset, prefix);
  } else if (type == "predict") {
    step =
        std::make_shared<steps::Predict>(*inputStep, parset, prefix, inputType);
  } else if (type == "grouppredict") {
    step = std::make_shared<steps::BdaGroupPredict>(*inputStep, parset, prefix);
  } else if (type == "idgpredict") {
    step = std::make_shared<steps::IDGPredict>(*inputStep, parset, prefix);
  } else if (type == "h5parmpredict") {
    step = std::make_shared<steps::H5ParmPredict>(inputStep, parset, prefix);
  } else if (type == "gaincal" || type == "calibrate") {
    step = std::make_shared<steps::GainCal>(*inputStep, parset, prefix);
  } else if (type == "upsample") {
    step = std::make_shared<steps::Upsample>(parset, prefix);
  } else if (type == "split" || type == "explode") {
    step = std::make_shared<steps::Split>(inputStep, parset, prefix);
  } else if (type == "ddecal") {
    if (inputType == Step::MsType::kRegular) {
      step = std::make_shared<steps::DDECal>(inputStep, parset, prefix);
    } else if (inputType == Step::MsType::kBda) {
      step = std::make_shared<steps::BdaDdeCal>(inputStep, parset, prefix);
    }
  } else if (type == "interpolate") {
    step = std::make_shared<steps::Interpolate>(inputStep, parset, prefix);
  } else if (type == "out" || type == "output" || type == "msout") {
    if (msName.empty())
      msName = casacore::Path(inputStep->msName()).absoluteName();
    step = makeOutputStep(inputStep, parset, prefix, msName, inputType);
  } else if (type == "python" || type == "pythondppp") {
    step = pythondp3::PyStep::create_instance(inputStep, parset, prefix);
  } else {
    // Maybe the step is defined in a dynamic library.
    step = findStepCtor(type)(inputStep, parset, prefix);
  }
  if (!step) {
    throw std::runtime_error("Could not create step of type '" + type + "'");
  }
  return step;
}

std::shared_ptr<Step> DP3::makeOutputStep(InputStep* reader,
                                          const common::ParameterSet& parset,
                                          const std::string& prefix,
                                          std::string& currentMSName,
                                          Step::MsType inputType) {
  std::shared_ptr<Step> step;
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
  switch (inputType) {
    case Step::MsType::kBda:
      if (doUpdate) {
        throw std::invalid_argument("No updater for BDA data.");
      } else {
        step = std::make_shared<MSBDAWriter>(reader, outName, parset, prefix);
        reader->setReadVisData(true);
      }
      break;
    case Step::MsType::kRegular:
      if (doUpdate) {
        // Create MSUpdater.
        // Take care the history is not written twice.
        // Note that if there is nothing to write, the updater won't do
        // anything.
        step = std::make_shared<MSUpdater>(reader, outName, parset, prefix,
                                           outName != currentMSName);
      } else {
        step = std::make_shared<MSWriter>(reader, outName, parset, prefix);
        reader->setReadVisData(true);
      }
      break;
  }
  casacore::Path pathOut(outName);
  currentMSName = pathOut.absoluteName();
  return step;
}

}  // namespace base
}  // namespace dp3
