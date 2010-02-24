//# NDPPP.cc: Program to execute steps like averaging and flagging on an MS
//# Copyright (C) 2010
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/DPStep.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <DPPP/MSReader.h>
#include <DPPP/MSWriter.h>
#include <DPPP/MSUpdater.h>
#include <DPPP/Averager.h>
#include <DPPP/MedFlagger.h>
#include <DPPP/ProgressMeter.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <iostream>
#include <stdexcept>
#include <libgen.h>

using namespace LOFAR::DPPP;
using namespace LOFAR;


DPStep::ShPtr makeSteps (const ParameterSet& parset)
{
  DPStep::ShPtr firstStep;
  DPStep::ShPtr lastStep;
  // Get input and output MS name.
  string inName  = parset.getString ("msin");
  string outName = parset.getString ("msout");
  // Get the steps.
  vector<string> steps = parset.getStringVector ("steps");
  // Currently the input MS must be given.
  // In the future it might be possible to have a simulation step instead.
  // Create MSReader step if input ms given.
  ASSERTSTR (!inName.empty(), "Name of input MS is not given");
  MSReader* reader = new MSReader (inName, parset, "msin.");
  firstStep = DPStep::ShPtr (reader);
  lastStep = firstStep;
  // Create the other steps.
  for (vector<string>::const_iterator iter = steps.begin();
       iter != steps.end(); ++iter) {
    string prefix(*iter + '.');
    string type = toLower(parset.getString (prefix+"type"));
    DPStep::ShPtr step;
    if (type == "average"  ||  type == "squash") {
      step = DPStep::ShPtr(new Averager (reader, parset, prefix));
    } else if (type == "madflagger") {
      step = DPStep::ShPtr(new MedFlagger (parset, prefix));
    } else {
      THROW (LOFAR::Exception, "DPPP step type " << type << " is unknown");
    }
    lastStep->setNextStep (step);
    lastStep = step;
    // Define as first step if not defined yet.
    if (!firstStep) {
      firstStep = step;
    }
  }
  // Find out how the data are averaged.
  AverageInfo avgInfo;
  DPStep::ShPtr step = firstStep;
  while (step) {
    step->updateAverageInfo (avgInfo);
    step = step->getNextStep();
  }
  // Create an updater step if an input MS was given; otherwise a writer.
  if (outName.empty()) {
    ASSERTSTR (avgInfo.nchanAvg() == 1  &&  avgInfo.ntimeAvg() == 1,
               "A new MS has to be given in msout if averaging is done");
    step = DPStep::ShPtr(new MSUpdater (reader, parset, ""));
  } else {
    step = DPStep::ShPtr(new MSWriter (reader, outName, avgInfo,
                                       parset, "msout."));
  }
  lastStep->setNextStep (step);
  lastStep = step;
  // Add a null step, so the last step can use getNextStep->process().
  DPStep::ShPtr nullStep(new NullStep());
  if (lastStep) {
    lastStep->setNextStep (nullStep);
  } else {
    firstStep = nullStep;
  }
  return firstStep;
}

int main(int argc, char *argv[])
{
  try
  {
    INIT_LOGGER(basename(argv[0]));
    // Create the steps.
    string parsetName("DPPP.parset");
    if (argc > 1) {
      parsetName = argv[1];
    }
    ParameterSet parset(parsetName);
    bool showProgress = parset.getBool ("showprogress", true);
    DPStep::ShPtr firstStep = makeSteps (parset);
    // Show the steps and determine AverageInfo again to get #times to process.
    AverageInfo avgInfo;
    DPStep::ShPtr step = firstStep;
    while (step) {
      step->updateAverageInfo (avgInfo);
      step->show (std::cout);
      step = step->getNextStep();
    }
    // Process until the end.
    uint ntodo = avgInfo.ntime() * avgInfo.ntimeAvg();
    std::cout << "Processing " << ntodo << " time slots ..." << std::endl;
    {
      ProgressMeter progress(0.0, ntodo, "DPPP", "Time slots processed",
                             "", "", true, 1);
      double ndone = 0;
      if (showProgress  &&  ntodo > 0) {
        progress.update (ndone, true);
      }
      DPBuffer buf;
      while (firstStep->process (buf)) {
        ++ndone;
        if (showProgress  &&  ntodo > 0) {
          progress.update (ndone, true);
        }
      }
    }
    // Finish the processing.
    std::cout << "Finishing processing ..." << std::endl;
    firstStep->finish();
    // Show the counts where needed.
    step = firstStep;
    while (step) {
      step->showCounts (std::cout);
      step = step->getNextStep();
    }
    // The destructors are called automatically at this point.
  } catch (std::exception& err) {
    std::cerr << "Error detected: " << err.what() << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "** PROBLEM **: Unhandled exception caught." << std::endl;
    return 2;
  }
}
