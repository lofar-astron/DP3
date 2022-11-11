// DPRun.cc: Class to run steps like averaging and flagging on an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/base/DP3.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "DPLogger.h"
#include "ProgressMeter.h"

#include <boost/algorithm/string.hpp>

#include "../steps/AntennaFlagger.h"
#include "../steps/AOFlaggerStep.h"
#include "../steps/ApplyBeam.h"
#include "../steps/ApplyCal.h"
#include "../steps/Averager.h"
#include "../steps/BDAAverager.h"
#include "../steps/BDAExpander.h"
#include "../steps/BdaGroupPredict.h"
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
#include "../steps/MsColumnReader.h"
#include "../steps/MSReader.h"
#include "../steps/MSUpdater.h"
#include "../steps/MSWriter.h"
#include "../steps/MultiMSReader.h"
#include "../steps/NullStep.h"
#include "../steps/NullStokes.h"
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

#include <dp3/common/Fields.h>
#include "../common/Timer.h"
#include "../common/StreamUtil.h"
#include "../steps/AntennaFlagger.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/OS/DynLib.h>
#include <casacore/casa/Utilities/Regex.h>

using dp3::steps::InputStep;
using dp3::steps::MSBDAWriter;
using dp3::steps::MSUpdater;
using dp3::steps::MSWriter;
using dp3::steps::OutputStep;
using dp3::steps::Step;
using dp3::common::operator<<;

namespace dp3 {
namespace base {

// Initialize the statics.
std::map<std::string, DP3::StepCreator*> DP3::theirStepMap;

/// Create an output step, either an MSWriter, MSUpdater or an MSBDAWriter
/// If no data are modified (for example if only count was done),
/// still an MSUpdater is created, but it will not write anything.
/// It reads the output name from the parset. If the prefix is "", it
/// reads msout or msout.name, otherwise it reads name from the output step
/// Create an updater step if an input MS was given; otherwise a writer.
/// Create an updater step only if needed (e.g. not if only count is done).
/// If the user specified an output MS name, a writer or updater is always
/// created If there is a writer, the reader needs to read the visibility
/// data.
static std::shared_ptr<OutputStep> MakeOutputStep(
    const common::ParameterSet& parset, const std::string& prefix,
    std::string& currentMSName, Step::MsType inputType) {
  std::shared_ptr<OutputStep> step;
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
  assert(!currentMSName.empty());
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
        step = std::make_shared<MSBDAWriter>(outName, parset, prefix);
      }
      break;
    case Step::MsType::kRegular:
      if (doUpdate) {
        // Create MSUpdater.
        // Take care the history is not written twice.
        // Note that if there is nothing to write, the updater won't do
        // anything.
        step = std::make_shared<MSUpdater>(outName, parset, prefix,
                                           outName != currentMSName);
      } else {
        step = std::make_shared<MSWriter>(outName, parset, prefix);
      }
      break;
  }
  casacore::Path pathOut(outName);
  currentMSName = pathOut.absoluteName();
  return step;
}

std::shared_ptr<Step> DP3::MakeSingleStep(const std::string& type,
                                          const common::ParameterSet& parset,
                                          const std::string& prefix,
                                          Step::MsType inputType) {
  std::shared_ptr<Step> step;
  if (type == "aoflagger" || type == "aoflag") {
    step = std::make_shared<steps::AOFlaggerStep>(parset, prefix);
  } else if (type == "averager" || type == "average" || type == "squash") {
    step = std::make_shared<steps::Averager>(parset, prefix);
  } else if (type == "bdaaverage" || type == "bdaaverager") {
    step = std::make_shared<steps::BDAAverager>(parset, prefix);
  } else if (type == "bdaexpander") {
    step = std::make_shared<steps::BDAExpander>(prefix);
  } else if (type == "madflagger" || type == "madflag") {
    step = std::make_shared<steps::MedFlagger>(parset, prefix);
  } else if (type == "preflagger" || type == "preflag") {
    step = std::make_shared<steps::PreFlagger>(parset, prefix);
  } else if (type == "antennaflagger") {
    step = std::make_shared<steps::AntennaFlagger>(parset, prefix);
  } else if (type == "uvwflagger" || type == "uvwflag") {
    step = std::make_shared<steps::UVWFlagger>(parset, prefix, inputType);
  } else if (type == "columnreader") {
    step = std::make_shared<steps::MsColumnReader>(parset, prefix);
  } else if (type == "counter" || type == "count") {
    step = std::make_shared<steps::Counter>(parset, prefix);
  } else if (type == "phaseshifter" || type == "phaseshift") {
    step = std::make_shared<steps::PhaseShift>(parset, prefix);
  } else if (type == "demixer" || type == "demix") {
    step = std::make_shared<steps::Demixer>(parset, prefix);
  } else if (type == "smartdemixer" || type == "smartdemix") {
    step = std::make_shared<steps::DemixerNew>(parset, prefix);
  } else if (type == "applybeam") {
    step = std::make_shared<steps::ApplyBeam>(parset, prefix);
  } else if (type == "stationadder" || type == "stationadd") {
    step = std::make_shared<steps::StationAdder>(parset, prefix);
  } else if (type == "scaledata") {
    step = std::make_shared<steps::ScaleData>(parset, prefix, inputType);
  } else if (type == "setbeam") {
    step = std::make_shared<steps::SetBeam>(parset, prefix);
  } else if (type == "filter") {
    step = std::make_shared<steps::Filter>(parset, prefix);
  } else if (type == "applycal" || type == "correct") {
    step = std::make_shared<steps::ApplyCal>(parset, prefix);
  } else if (type == "nullstokes") {
    step = std::make_shared<steps::NullStokes>(parset, prefix);
  } else if (type == "predict") {
    step = std::make_shared<steps::Predict>(parset, prefix, inputType);
  } else if (type == "idgpredict") {
    step = std::make_shared<steps::IDGPredict>(parset, prefix);
  } else if (type == "upsample") {
    step = std::make_shared<steps::Upsample>(parset, prefix);
  } else if (type == "interpolate") {
    step = std::make_shared<steps::Interpolate>(parset, prefix);
  } else if (type == "null") {
    step = std::make_shared<steps::NullStep>();
  } else if (type == "smartdemixer" || type == "smartdemix") {
    step = std::make_shared<steps::DemixerNew>(parset, prefix);
  } else if (type == "grouppredict") {
    step = std::make_shared<steps::BdaGroupPredict>(parset, prefix);
  } else if (type == "idgpredict") {
    step = std::make_shared<steps::IDGPredict>(parset, prefix);
  } else if (type == "h5parmpredict") {
    step = std::make_shared<steps::H5ParmPredict>(parset, prefix);
  } else if (type == "gaincal" || type == "calibrate") {
    step = std::make_shared<steps::GainCal>(parset, prefix);
  } else if (type == "demixer" || type == "demix") {
    step = std::make_shared<steps::Demixer>(parset, prefix);
  } else if (type == "python" || type == "pythondppp") {
    step = pythondp3::PyStep::create_instance(parset, prefix);
  } else if (type == "null") {
    step = std::make_shared<steps::NullStep>();
  }
  return step;
}

static std::shared_ptr<Step> MakeSingleStep(const std::string& type,
                                            InputStep* inputStep,
                                            const common::ParameterSet& parset,
                                            const std::string& prefix,
                                            std::string& msName,
                                            Step::MsType inputType) {
  std::shared_ptr<Step> step =
      DP3::MakeSingleStep(type, parset, prefix, inputType);
  if (step) {
    return step;
  }
  if (type == "split" || type == "explode") {
    step = std::make_shared<steps::Split>(inputStep, parset, prefix);
  } else if (type == "ddecal") {
    if (inputType == Step::MsType::kRegular) {
      step = std::make_shared<steps::DDECal>(inputStep, parset, prefix);
    } else if (inputType == Step::MsType::kBda) {
      step = std::make_shared<steps::BdaDdeCal>(parset, prefix);
    }
  } else if (type == "out" || type == "output" || type == "msout") {
    step = MakeOutputStep(parset, prefix, msName, inputType);
  } else {
    // Maybe the step is defined in a dynamic library.
    step = DP3::FindStepCreator(type)(parset, prefix);
  }
  if (!step) {
    throw std::runtime_error("Could not create step of type '" + type + "'");
  }
  return step;
}

void DP3::RegisterStepCreator(const std::string& type, StepCreator* func) {
  theirStepMap[type] = func;
}

DP3::StepCreator* DP3::FindStepCreator(const std::string& type) {
  std::map<std::string, StepCreator*>::const_iterator iter =
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
  casacore::DynLib dl(libname, "libDP3_", "register_" + libname, false);
  if (dl.getHandle()) {
    // See if registered now.
    iter = theirStepMap.find(type);
    if (iter != theirStepMap.end()) {
      return iter->second;
    }
  }

  throw std::runtime_error(
      "Step type " + type + " is unknown and no shared library lib" + libname +
      " or libDP3_" + libname + " found in (DY)LD_LIBRARY_PATH");
}

dp3::common::Fields DP3::GetChainRequiredFields(
    std::shared_ptr<Step> first_step) {
  std::shared_ptr<Step> last_step;
  for (std::shared_ptr<Step> step = first_step; step;
       step = step->getNextStep()) {
    last_step = step;
  }

  dp3::common::Fields overall_required_fields;
  for (Step* step = last_step.get(); step; step = step->getPrevStep()) {
    overall_required_fields.UpdateRequirements(step->getRequiredFields(),
                                               step->getProvidedFields());
    if (step == first_step.get()) break;
  }

  return overall_required_fields;
}

dp3::common::Fields DP3::SetChainProvidedFields(
    std::shared_ptr<Step> first_step, common::Fields provided_fields) {
  for (std::shared_ptr<Step> step = first_step; step;
       step = step->getNextStep()) {
    OutputStep* output_step = dynamic_cast<OutputStep*>(step.get());
    if (output_step) {
      output_step->SetFieldsToWrite(provided_fields);
      provided_fields = common::Fields();
    } else {
      provided_fields |= step->getProvidedFields();
    }
  }
  return provided_fields;
}

void DP3::Execute(const string& parsetName, int argc, char* argv[]) {
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
  std::shared_ptr<InputStep> firstStep = MakeMainSteps(parset);

  Step::ShPtr step = firstStep;
  Step::ShPtr lastStep;
  while (step) {
    step = step->getNextStep();
  }

  // Tell the reader which fields must be read.
  firstStep->setFieldsToRead(GetChainRequiredFields(firstStep));

  // Call updateInfo()
  DPInfo dpInfo;
  if (numThreads > 0) {
    dpInfo.setNThreads(numThreads);
  }
  dpInfo = firstStep->setInfo(dpInfo);

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
      if (checkparset != 0)
        throw std::runtime_error("Unused parset keywords found");
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
  timer.show(ostr, "Total DP3 time");
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

std::shared_ptr<InputStep> DP3::MakeMainSteps(
    const common::ParameterSet& parset) {
  std::shared_ptr<InputStep> input_step = InputStep::CreateReader(parset);
  std::shared_ptr<Step> last_step = input_step;

  std::shared_ptr<Step> step = MakeStepsFromParset(
      parset, "", "steps", *input_step, false, input_step->outputs());
  if (step) {
    input_step->setNextStep(step);

    while (step->getNextStep()) step = step->getNextStep();
    last_step = step;
  }

  // Determine the provided fields of the series of steps. When provided
  // fields is non-empty, create an output step that writes those fields.
  const common::Fields provided_fields = SetChainProvidedFields(input_step);

  // Check if the last step is an output step.
  const bool ends_with_output_step = dynamic_cast<OutputStep*>(last_step.get());

  if (!ends_with_output_step) {
    // Check if an output step is needed because of the parset.
    const std::string msOutName = parset.getString(
        parset.isDefined("msout.name") ? "msout.name" : "msout");

    if (!msOutName.empty() || provided_fields != common::Fields()) {
      std::string msName = casacore::Path(input_step->msName()).absoluteName();
      std::shared_ptr<OutputStep> output_step =
          MakeOutputStep(parset, "msout.", msName, last_step->outputs());
      output_step->SetFieldsToWrite(provided_fields);
      last_step->setNextStep(output_step);
      last_step = output_step;
    }
  }

  // Add a null step, so the last step can use getNextStep->process().
  // Split may not have a next step (Split::setNextStep throws).
  if (!dynamic_cast<steps::Split*>(last_step.get())) {
    last_step->setNextStep(std::make_shared<steps::NullStep>());
  }
  return input_step;
}

Step::ShPtr DP3::MakeStepsFromParset(const common::ParameterSet& parset,
                                     const std::string& prefix,
                                     const std::string& step_names_key,
                                     InputStep& inputStep, bool terminateChain,
                                     Step::MsType initial_step_output) {
  std::string msName = casacore::Path(inputStep.msName()).absoluteName();
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
    Step::ShPtr step = dp3::base::MakeSingleStep(type, &inputStep, parset,
                                                 prefix, msName, inputType);

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

  if (terminateChain && lastStep) {
    // Add a null step, so the last step can use getNextStep->process().
    lastStep->setNextStep(std::make_shared<steps::NullStep>());
  }

  return firstStep;
}

}  // namespace base
}  // namespace dp3
