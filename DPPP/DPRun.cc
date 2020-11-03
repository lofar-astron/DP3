// DPRun.cc: Class to run steps like averaging and flagging on an MS
// Copyright (C) 2010
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$
//
// @author Ger van Diepen

#include "DPRun.h"

#include <boost/algorithm/string.hpp>

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "Averager.h"
#include "BDAAverager.h"
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
  DPStep::ShPtr firstStep = makeSteps(parset, "", 0);

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
  firstStep->setInfo(std::move(dpInfo));

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

DPStep::ShPtr DPRun::makeSteps(const ParameterSet& parset, const string& prefix,
                               DPInput* reader) {
  DPStep::ShPtr firstStep;
  DPStep::ShPtr lastStep;
  if (!reader) {
    std::unique_ptr<DPInput> new_reader = DPInput::CreateReader(parset, prefix);
    reader = new_reader.get();
    firstStep = std::move(new_reader);
  }

  casacore::Path pathIn(reader->msName());
  casacore::String currentMSName(pathIn.absoluteName());

  // Create the other steps.
  std::vector<string> steps = parset.getStringVector(prefix + "steps");
  lastStep = firstStep;
  DPStep::ShPtr step;
  bool needsOutputStep = true;

  for (std::vector<string>::const_iterator iter = steps.begin();
       iter != steps.end(); ++iter) {
    string prefix(*iter + '.');
    // The alphabetic part of the name is the default step type.
    // This allows names like average1, out3.
    string defaulttype = (*iter);
    while (defaulttype.size() > 0 && std::isdigit(*defaulttype.rbegin())) {
      defaulttype.resize(defaulttype.size() - 1);
    }

    // If no explicit output step is given as last step, one will be added
    // with the msout. prefix
    needsOutputStep = true;

    string type = parset.getString(prefix + "type", defaulttype);
    boost::algorithm::to_lower(type);
    // Define correct name for AOFlagger synonyms.
    if (type == "aoflagger" || type == "aoflag") {
      step = std::make_shared<AOFlaggerStep>(reader, parset, prefix);
    } else if (type == "averager" || type == "average" || type == "squash") {
      step = std::make_shared<Averager>(reader, parset, prefix);
    } else if (type == "bdaaverager") {
      step = std::make_shared<BDAAverager>(*reader, parset, prefix);
    } else if (type == "madflagger" || type == "madflag") {
      step = std::make_shared<MedFlagger>(reader, parset, prefix);
    } else if (type == "preflagger" || type == "preflag") {
      step = std::make_shared<PreFlagger>(reader, parset, prefix);
    } else if (type == "uvwflagger" || type == "uvwflag") {
      step = std::make_shared<UVWFlagger>(reader, parset, prefix);
    } else if (type == "counter" || type == "count") {
      step = std::make_shared<Counter>(reader, parset, prefix);
    } else if (type == "phaseshifter" || type == "phaseshift") {
      step = std::make_shared<PhaseShift>(reader, parset, prefix);
    } else if (type == "demixer" || type == "demix") {
      step = std::make_shared<Demixer>(reader, parset, prefix);
    } else if (type == "smartdemixer" || type == "smartdemix") {
      step = std::make_shared<DemixerNew>(reader, parset, prefix);
    } else if (type == "applybeam") {
      step = std::make_shared<ApplyBeam>(reader, parset, prefix);
    } else if (type == "stationadder" || type == "stationadd") {
      step = std::make_shared<StationAdder>(reader, parset, prefix);
    } else if (type == "scaledata") {
      step = std::make_shared<ScaleData>(reader, parset, prefix);
    } else if (type == "setbeam") {
      step = std::make_shared<SetBeam>(reader, parset, prefix);
    } else if (type == "filter") {
      step = std::make_shared<Filter>(reader, parset, prefix);
    } else if (type == "applycal" || type == "correct") {
      step = std::make_shared<ApplyCal>(reader, parset, prefix);
    } else if (type == "predict") {
      step = std::make_shared<Predict>(reader, parset, prefix);
    } else if (type == "idgpredict") {
      step = std::make_shared<IDGPredict>(*reader, parset, prefix);
    } else if (type == "h5parmpredict") {
      step = std::make_shared<H5ParmPredict>(reader, parset, prefix);
    } else if (type == "gaincal" || type == "calibrate") {
      step = std::make_shared<GainCal>(reader, parset, prefix);
    } else if (type == "upsample") {
      step = std::make_shared<Upsample>(reader, parset, prefix);
    } else if (type == "split" || type == "explode") {
      step = std::make_shared<Split>(reader, parset, prefix);
      needsOutputStep = false;
    } else if (type == "ddecal") {
      step = std::make_shared<DDECal>(reader, parset, prefix);
    } else if (type == "interpolate") {
      step = std::make_shared<Interpolate>(reader, parset, prefix);
    } else if (type == "out" || type == "output" || type == "msout") {
      step = makeOutputStep(reader, parset, prefix, currentMSName,
                            lastStep->outputs() == DPStep::MSType::BDA);
      needsOutputStep = false;
    } else {
      // Maybe the step is defined in a dynamic library.
      step = findStepCtor(type)(reader, parset, prefix);
    }

    if (!step->accepts(lastStep->outputs())) {
      throw std::invalid_argument("Step " + type +
                                  " is incompatible with input data.");
    }

    if (lastStep) {
      lastStep->setNextStep(step);
    }
    lastStep = step;
    // Define as first step if not defined yet.
    if (!firstStep) {
      firstStep = step;
    }
  }
  // Add an output step if not explicitly added in steps (unless last step is a
  // 'split' step)
  if (needsOutputStep) {
    step = makeOutputStep(reader, parset, "msout.", currentMSName,
                          lastStep->outputs() == DPStep::MSType::BDA);
    lastStep->setNextStep(step);
    lastStep = step;
  }

  // Add a null step, so the last step can use getNextStep->process().
  auto nullStep = std::make_shared<NullStep>();
  if (lastStep) {
    lastStep->setNextStep(nullStep);
  } else {
    firstStep = nullStep;
  }
  return firstStep;
}

DPStep::ShPtr DPRun::makeOutputStep(DPInput* reader, const ParameterSet& parset,
                                    const string& prefix,
                                    casacore::String& currentMSName,
                                    const bool& isBDA) {
  DPStep::ShPtr step;
  casacore::String outName;
  bool doUpdate = false;
  if (prefix == "msout.") {
    // The last output step.
    outName = parset.getString("msout.name", "");
    if (outName.empty()) {
      outName = parset.getString("msout");
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
    if (currentMSName == pathOut.absoluteName()) {
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
