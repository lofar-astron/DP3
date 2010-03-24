//# DPRun.cc: Class to run steps like averaging and flagging on an MS
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
#include <DPPP/DPRun.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <DPPP/MSReader.h>
#include <DPPP/MSWriter.h>
#include <DPPP/MSUpdater.h>
#include <DPPP/Averager.h>
#include <DPPP/MedFlagger.h>
#include <DPPP/PreFlagger.h>
#include <DPPP/UVWFlagger.h>
#include <DPPP/Counter.h>
#include <DPPP/ProgressMeter.h>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>
#include <casa/OS/Timer.h>

namespace LOFAR {
  namespace DPPP {

    void DPRun::execute (const string& parsetName)
    {
      casa::Timer timer;
      NSTimer nstimer;
      nstimer.start();
      ParameterSet parset(parsetName);
      bool showProgress = parset.getBool ("showprogress", true);
      bool showTimings  = parset.getBool ("showtimings", true);
      DPStep::ShPtr firstStep = makeSteps (parset);
      // Show the steps and determine AverageInfo again to get #times
      // to process.
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
        ProgressMeter* progress = 0;
        if (showProgress) {
          progress = new ProgressMeter(0.0, ntodo, "NDPPP",
                                       "Time slots processed",
                                       "", "", true, 1);
        }
        double ndone = 0;
        if (showProgress  &&  ntodo > 0) {
          progress->update (ndone, true);
        }
        DPBuffer buf;
        while (firstStep->process (buf)) {
          ++ndone;
          if (showProgress  &&  ntodo > 0) {
            progress->update (ndone, true);
          }
        }
        delete progress;
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
      // Show the overall timer.
      nstimer.stop();
      double duration = nstimer.getElapsed();
      cout << endl;
      timer.show ("Total NDPPP time");
      if (showTimings) {
        // Show the timings per step.
        step = firstStep;
        while (step) {
          step->showTimings (std::cout, duration);
          step = step->getNextStep();
        }
      }
      // The destructors are called automatically at this point.
    }

    DPStep::ShPtr DPRun::makeSteps (const ParameterSet& parset)
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
        if (type == "averager"  ||  type == "average"  ||  type == "squash") {
          step = DPStep::ShPtr(new Averager (reader, parset, prefix));
        } else if (type == "madflagger"  ||  type == "madflag") {
          step = DPStep::ShPtr(new MedFlagger (reader, parset, prefix));
        } else if (type == "preflagger"  ||  type == "preflag") {
          step = DPStep::ShPtr(new PreFlagger (reader, parset, prefix));
        } else if (type == "uvwflagger"  ||  type == "uvwflag") {
          step = DPStep::ShPtr(new UVWFlagger (reader, parset, prefix));
        } else if (type == "counter"  ||  type == "count") {
          step = DPStep::ShPtr(new Counter (reader, parset, prefix));
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
        step = DPStep::ShPtr(new MSUpdater (reader, parset, "msout."));
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

  } //# end namespace
}
