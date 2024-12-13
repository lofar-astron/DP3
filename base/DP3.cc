// DPRun.cc: Class to run steps like averaging and flagging on an MS
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/base/DP3.h>

#include <aocommon/checkblas.h>
#include <aocommon/logger.h>
#include <aocommon/threadpool.h>

#include <boost/algorithm/string.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "ProgressMeter.h"
#include "SkyModelCache.h"
#include "Version.h"

#include "../steps/AntennaFlagger.h"
#include "../steps/AOFlaggerStep.h"
#include "../steps/ApplyBeam.h"
#include "../steps/ApplyCal.h"
#include "../steps/Averager.h"
#include "../steps/BDAAverager.h"
#include "../steps/BdaExpander.h"
#include "../steps/BdaGroupPredict.h"
#include "../steps/Clipper.h"
#include "../steps/Counter.h"
#include "../steps/DDECal.h"
#include "../steps/BdaDdeCal.h"
#include "../steps/Demixer.h"
#include "../steps/Filter.h"
#include "../steps/GainCal.h"
#include "../steps/H5ParmPredict.h"
#include "../steps/IDGPredict.h"
#include "../steps/IDGImager.h"
#include "../steps/Interpolate.h"
#include "../steps/MadFlagger.h"
#include "../steps/MSBDAWriter.h"
#include "../steps/MsColumnReader.h"
#include "../steps/MSUpdater.h"
#include "../steps/MSWriter.h"
#include "../steps/WSCleanWriter.h"
#include "../steps/NullStep.h"
#include "../steps/NullStokes.h"
#include "../steps/PhaseShift.h"
#include "../steps/Predict.h"
#include "../steps/PreFlagger.h"
#include "../steps/SagecalPredict.h"
#include "../steps/ScaleData.h"
#include "../steps/SetBeam.h"
#include "../steps/Split.h"
#include "../steps/StationAdder.h"
#include "../steps/UVWFlagger.h"
#include "../steps/Upsample.h"
#include "../steps/WGridderPredict.h"
#include "../steps/FlagTransfer.h"

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
  if (outName.empty() || outName == ".") {
    // currentMSName is empty when creating sub-steps, e.g, in DDECal and Split.
    if (currentMSName.empty()) {
      throw std::runtime_error(
          "In a series of steps that are part of another step, the first "
          "output step must have a measurement set name.");
    }
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

std::shared_ptr<Step> MakeSingleStep(const std::string& type,
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
    step = std::make_shared<steps::BdaExpander>(prefix);
  } else if (type == "madflagger" || type == "madflag") {
    step = std::make_shared<steps::MadFlagger>(parset, prefix);
  } else if (type == "preflagger" || type == "preflag") {
    step = std::make_shared<steps::PreFlagger>(parset, prefix);
  } else if (type == "antennaflagger" || type == "antflag") {
    step = std::make_shared<steps::AntennaFlagger>(parset, prefix);
  } else if (type == "uvwflagger" || type == "uvwflag") {
    step = std::make_shared<steps::UVWFlagger>(parset, prefix, inputType);
  } else if (type == "clipper") {
    step = std::make_shared<steps::Clipper>(parset, prefix);
  } else if (type == "columnreader") {
    step = std::make_shared<steps::MsColumnReader>(parset, prefix);
  } else if (type == "counter" || type == "count") {
    step = std::make_shared<steps::Counter>(parset, prefix);
  } else if (type == "phaseshifter" || type == "phaseshift") {
    step = std::make_shared<steps::PhaseShift>(parset, prefix);
  } else if (type == "demixer" || type == "demix") {
    step = std::make_shared<steps::Demixer>(parset, prefix);
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
  } else if (type == "wgridderpredict") {
    step = std::make_shared<steps::WGridderPredict>(parset, prefix);
  } else if (type == "idgpredict") {
    step = std::make_shared<steps::IDGPredict>(parset, prefix);
  } else if (type == "idgimager") {
    step = std::make_shared<steps::IDGImager>(parset, prefix);
  } else if (type == "upsample") {
    step = std::make_shared<steps::Upsample>(parset, prefix);
  } else if (type == "interpolate") {
    step = std::make_shared<steps::Interpolate>(parset, prefix);
  } else if (type == "grouppredict") {
    step = std::make_shared<steps::BdaGroupPredict>(parset, prefix);
  } else if (type == "sagecalpredict") {
    step = std::make_shared<steps::SagecalPredict>(parset, prefix);
  } else if (type == "h5parmpredict") {
    step = std::make_shared<steps::H5ParmPredict>(parset, prefix);
  } else if (type == "gaincal" || type == "calibrate") {
    step = std::make_shared<steps::GainCal>(parset, prefix);
  } else if (type == "python" || type == "pythondppp") {
    step = pythondp3::PyStep::create_instance(parset, prefix);
  } else if (type == "split" || type == "explode") {
    step = std::make_shared<steps::Split>(parset, prefix);
  } else if (type == "ddecal") {
    if (inputType == Step::MsType::kRegular) {
      step = std::make_shared<steps::DDECal>(parset, prefix);
    } else if (inputType == Step::MsType::kBda) {
      step = std::make_shared<steps::BdaDdeCal>(parset, prefix);
    }
  } else if (type == "null") {
    step = std::make_shared<steps::NullStep>();
  } else if (type == "wscleanwriter") {
    step = std::make_shared<steps::WSCleanWriter>(parset, prefix);
  } else if (type == "flagtransfer") {
    step = std::make_shared<steps::FlagTransfer>(parset, prefix);
  }
  return step;
}

dp3::common::Fields GetChainRequiredFields(std::shared_ptr<Step> first_step) {
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

dp3::common::Fields SetChainProvidedFields(std::shared_ptr<Step> first_step,
                                           common::Fields provided_fields) {
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

void Execute(const std::string& parsetName, int argc, char* argv[]) {
  casacore::Timer timer;
  common::NSTimer nstimer;
  nstimer.start();
  common::ParameterSet parset;
  if (!parsetName.empty()) {
    parset.adoptFile(parsetName);
  }
  // Adopt possible parameters given at the command line.
  parset.adoptArgv(argc, argv);  ///< works fine if argc==0 and argv==0

  // Immediately initialize logger such that output will follow requested
  // verbosity
  if (parset.isDefined("verbosity")) {
    aocommon::Logger::SetVerbosity(aocommon::StringToLogVerbosityLevel(
        parset.getString("verbosity", "normal")));
  }
  aocommon::Logger::SetLogTime(parset.getBool("time_logging", false));
  aocommon::Logger::SetLogMemory(parset.getBool("memory_logging", false));

  bool showProgress = parset.getBool("showprogress", true);
  bool showTimings = parset.getBool("showtimings", true);
  // checkparset is an integer parameter now, but accepts a bool as well
  // for backward compatibility.
  int checkparset = 0;
  try {
    checkparset = parset.getInt("checkparset", 0);
  } catch (...) {
    aocommon::Logger::Warn
        << "Parameter checkparset should be an integer value\n";
    checkparset = parset.getBool("checkparset") ? 1 : 0;
  }

  bool showcounts = parset.getBool("showcounts", true);

  size_t n_threads = parset.getInt("numthreads", 0);
  if (n_threads == 0) n_threads = aocommon::system::ProcessorCount();
  aocommon::ThreadPool::GetInstance().SetNThreads(n_threads);
  Step::SetThreadingIsInitialized();
  aocommon::Logger::Debug << "DP3 started with " << n_threads << " threads.\n";

  // Create the steps, link them together
  std::shared_ptr<InputStep> firstStep = MakeMainSteps(parset);

  // Call updateInfo() on all steps
  DPInfo dpInfo = firstStep->setInfo(DPInfo());

  // Show the steps.
  std::shared_ptr<Step> step = firstStep;
  std::shared_ptr<Step> lastStep;
  while (step) {
    std::ostringstream os;
    step->show(os);
    aocommon::Logger::Info << os.str();
    lastStep = step;
    step = step->getNextStep();
  }
  if (checkparset >= 0) {
    // Show unused parameters (might be misspelled).
    std::vector<std::string> unused = parset.unusedKeys();
    if (!unused.empty()) {
      aocommon::Logger::Warn
          << "\n*** WARNING: the following parset keywords were not used ***"
          << "\n             maybe they are misspelled"
          << "\n";
      for (const std::string& s : unused)
        aocommon::Logger::Warn << "    - " << s << '\n';
      aocommon::Logger::Warn << '\n';
      if (checkparset != 0)
        throw std::runtime_error("Unused parset keywords found");
    }
  }

  // All steps should have finished reading the sky model, so clear the cache
  SkyModelCache::GetInstance().Clear();

  // Process until the end.
  unsigned int ntodo = firstStep->getInfo().ntime();
  aocommon::Logger::Info << "Processing " << ntodo << " time slots ...\n";
  if (showProgress) {
    double ndone = 0;
    ProgressMeter progress(ndone, ntodo, "DP3", "Time slots processed", "", "",
                           true, 1);
    if (ntodo > 0) progress.update(ndone, true);
    while (firstStep->process(std::make_unique<DPBuffer>())) {
      ++ndone;
      if (ntodo > 0) progress.update(ndone, true);
    }
  } else {
    while (firstStep->process(std::make_unique<DPBuffer>())) {
      // do nothing
    }
  }
  // Finish the processing.
  aocommon::Logger::Info << "Finishing processing ...\n";
  firstStep->finish();

  // Show the counts where needed.
  if (showcounts) {
    step = firstStep;
    while (step) {
      std::ostringstream os;
      step->showCounts(os);
      aocommon::Logger::Info << os.str();
      step = step->getNextStep();
    }
  }
  // Show the overall timer.
  nstimer.stop();
  double duration = nstimer.getElapsed();
  std::ostringstream ostr;
  ostr << '\n';
  // Output special line for pipeline use.
  aocommon::Logger::Debug << "Start timer output\n";
  timer.show(ostr, "Total DP3 time");
  aocommon::Logger::Info << ostr.str();
  if (showTimings) {
    // Show the timings per step.
    step = firstStep;
    while (step) {
      std::ostringstream os;
      step->showTimings(os, duration);
      if (!os.str().empty()) {
        aocommon::Logger::Info << os.str();
      }
      step = step->getNextStep();
    }
  }
  aocommon::Logger::Debug << "End timer output\n";
  // The destructors are called automatically at this point.
}

void ShowUsage() {
  aocommon::Logger::Info
      << "Usage: DP3 [-v] [parsetfile] [parsetkeys...]\n"
         "  parsetfile: a file containing one parset key=value pair per line\n"
         "  parsetkeys: any number of parset key=value pairs, e.g. "
         "msin=my.MS\n\n"
         "If both a file and command-line keys are specified, the keys on the "
         "command\n"
         "line override those in the file.\n"
         "If no arguments are specified, the program tries to read "
         "\"DP3.parset\",\n"
         "\"NDPPP.parset\" or \"DPPP.parset\" as a default.\n"
         "-v will show version info and exit.\n"
         "Documentation is at: https://dp3.readthedocs.io\n";
}

void ExecuteFromCommandLine(int argc, char* argv[]) {
  check_openblas_multithreading();

  // Get the name of the parset file.
  if (argc > 1) {
    string param = argv[1];
    if (param == "--help" || param == "-help" || param == "-h" ||
        param == "--usage" || param == "-usage") {
      ShowUsage();
      return;
    } else if (param == "-v" || param == "--version") {
      aocommon::Logger::Info << DP3Version::AsString(true) << '\n';
      return;
    }
  }

  std::string parsetName;
  if (argc > 1 && std::string(argv[1]).find('=') == std::string::npos) {
    // First argument is parset name (except if it's a key-value pair)
    parsetName = argv[1];
  } else if (argc == 1) {
    // No arguments given: try to load [N]DPPP.parset
    if (std::filesystem::exists("DP3.parset")) {
      parsetName = "DP3.parset";
    } else if (std::filesystem::exists("DPPP.parset")) {
      parsetName = "DPPP.parset";
    } else if (std::filesystem::exists("NDPPP.parset")) {
      parsetName = "NDPPP.parset";
    } else {  // No default file, show usage and exit
      ShowUsage();
      return;
    }
  }

  // Execute the parset file.
  Execute(parsetName, argc, argv);
}

std::shared_ptr<InputStep> MakeMainSteps(const common::ParameterSet& parset) {
  std::shared_ptr<InputStep> input_step = InputStep::CreateReader(parset);
  std::shared_ptr<Step> last_step = input_step;

  // Create the second and later steps, as requested by the parset. The chain
  // is not terminated by a null step yet.
  const std::string ms_name =
      casacore::Path(input_step->msName()).absoluteName();
  std::shared_ptr<Step> step = MakeStepsFromParset(
      parset, "", "steps", ms_name, false, input_step->outputs());
  if (step) {
    input_step->setNextStep(step);

    while (step->getNextStep()) step = step->getNextStep();
    last_step = step;
  }

  // Determine the provided fields of the series of steps. When provided_fields
  // is non-empty, create an output step that writes those fields.
  const common::Fields provided_fields = SetChainProvidedFields(input_step);

  // Check if the last step is an output step. If not, add it when necessary
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

  // Tell the reader which fields must be read.
  input_step->setFieldsToRead(
      GetChainRequiredFields(input_step->getNextStep()));

  return input_step;
}

std::shared_ptr<Step> MakeStepsFromParset(const common::ParameterSet& parset,
                                          const std::string& prefix,
                                          const std::string& step_names_key,
                                          const std::string& input_ms_name,
                                          bool terminateChain,
                                          Step::MsType initial_step_output) {
  std::string msName = input_ms_name;
  const std::vector<std::string> stepNames =
      parset.getStringVector(prefix + step_names_key);

  std::shared_ptr<Step> firstStep;
  std::shared_ptr<Step> lastStep;
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
    std::shared_ptr<Step> step =
        MakeSingleStep(type, parset, prefix, inputType);
    if (!step && (type == "out" || type == "output" || type == "msout")) {
      step = MakeOutputStep(parset, prefix, msName, inputType);
    }
    if (!step) {
      throw std::runtime_error("Could not create step of type '" + type + "'");
    }

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

size_t GetNThreads() { return aocommon::ThreadPool::GetInstance().NThreads(); }

void SetNThreads(size_t n_threads) {
  aocommon::ThreadPool::GetInstance().SetNThreads(n_threads);
  Step::SetThreadingIsInitialized();
}

}  // namespace base
}  // namespace dp3
