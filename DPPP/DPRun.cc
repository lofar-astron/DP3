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

#include "DPRun.h"

#include <boost/algorithm/string.hpp>

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "Averager.h"
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
#include "MSReader.h"
#include "MSUpdater.h"
#include "MSWriter.h"
#include "MultiMSReader.h"
#include "PhaseShift.h"
#include "Predict.h"
#include "PreFlagger.h"
#include "ProgressMeter.h"
#include "ScaleData.h"
#include "Split.h"
#include "StationAdder.h"
#include "UVWFlagger.h"
#include "Upsample.h"

#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"

#include "../AOFlaggerStep/AOFlaggerStep.h"

#include "../DDECal/DDECal.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/OS/DynLib.h>
#include <casacore/casa/Utilities/Regex.h>

namespace DP3 {
  namespace DPPP {

    // Initialize the statics.
    std::map<std::string, DPRun::StepCtor*> DPRun::theirStepMap;

    void DPRun::registerStepCtor (const std::string& type, StepCtor* func)
    {
      theirStepMap[type] = func;
    }

    DPRun::StepCtor* DPRun::findStepCtor (const std::string& type)
    {
      std::map<std::string,StepCtor*>::const_iterator iter =
        theirStepMap.find (type);
      if (iter != theirStepMap.end()) {
        return iter->second;
      }
      // Try to load the step from a dynamic library with that name
      // (in lowercase).
      // A dot can be used to have a specific library name (so multiple
      // steps can use the same shared library).
      std::string libname(type);
      boost::algorithm::to_lower(libname);
      string::size_type pos = libname.find_first_of (".");
      if (pos != string::npos) {
        libname = libname.substr (0, pos);
      }
      // Try to load and initialize the dynamic library.
      casacore::DynLib dl(libname, string("libdppp_"), "register_"+libname, false);
      if (dl.getHandle()) {
        // See if registered now.
        iter = theirStepMap.find (type);
        if (iter != theirStepMap.end()) {
          return iter->second;
        }
      }
      throw Exception("Step type " + type +
            " is unknown and no shared library lib" + libname + " or libdppp_" +
            libname + " found in (DY)LD_LIBRARY_PATH");
    }


    void DPRun::execute (const string& parsetName, int argc, char* argv[])
    {
      casacore::Timer timer;
      NSTimer nstimer;
      nstimer.start();
      ParameterSet parset;
      if (! parsetName.empty()) {
        parset.adoptFile (parsetName);
      }
      // Adopt possible parameters given at the command line.
      parset.adoptArgv (argc, argv); //# works fine if argc==0 and argv==0
      DPLogger::useLogger = parset.getBool ("uselogger", false);
      bool showProgress   = parset.getBool ("showprogress", true);
      bool showTimings    = parset.getBool ("showtimings", true);
      // checkparset is an integer parameter now, but accept a bool as well
      // for backward compatibility.
      int checkparset = 0;
      try {
        checkparset = parset.getInt ("checkparset", 0);
      } catch (...) {
        DPLOG_WARN_STR ("Parameter checkparset should be an integer value");
        checkparset = parset.getBool ("checkparset") ? 1:0;
      }

      bool showcounts = parset.getBool ("showcounts", true);

      uint numThreads = parset.getInt("numthreads", ThreadPool::NCPUs());

      // Create the steps, link them together
      DPStep::ShPtr firstStep = makeSteps (parset, "", 0);

      DPStep::ShPtr step = firstStep;
      DPStep::ShPtr lastStep;
      while (step) {
        step = step->getNextStep();
      }
      
      // Call updateInfo()
      DPInfo dpInfo;
      dpInfo.setNThreads(numThreads);
      firstStep->setInfo (std::move(dpInfo));

      // Show the steps.  
      step = firstStep;
      while (step) {
        std::ostringstream os;
        step->show (os);
        DPLOG_INFO (os.str(), true);
        lastStep = step;
        step = step->getNextStep();
      }
      if (checkparset >= 0) {
        // Show unused parameters (might be misspelled).
        std::vector<std::string> unused = parset.unusedKeys();
        if (! unused.empty()) {
          DPLOG_WARN_STR
            (
                "\n*** WARNING: the following parset keywords were not used ***"
             << "\n             maybe they are misspelled"
             << "\n    " << unused << std::endl);
          if (checkparset!=0)
            throw Exception("Unused parset keywords found");
        }
      }
      // Process until the end.
      uint ntodo = firstStep->getInfo().ntime();
      DPLOG_INFO_STR ("Processing " << ntodo << " time slots ...");
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
      DPLOG_INFO_STR ("Finishing processing ...");
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
          step->showCounts (os);
          DPLOG_INFO (os.str(), true);
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
      timer.show (ostr, "Total NDPPP time");
      DPLOG_INFO (ostr.str(), true);
      if (showTimings) {
        // Show the timings per step.
        step = firstStep;
        while (step) {
          std::ostringstream os;
          step->showTimings (os, duration);
        if (! os.str().empty()) {
          DPLOG_INFO (os.str(), true);
        }
          step = step->getNextStep();
        }
      }
      if (DPLogger::useLogger) {
        ostr << "End timer output\n";
      }
      // The destructors are called automatically at this point.
    }

    DPStep::ShPtr DPRun::makeSteps (const ParameterSet& parset,
                                    const string& prefix,
                                    DPInput* reader)
    {
      DPStep::ShPtr firstStep;
      DPStep::ShPtr lastStep;
      if (!reader) {
        // Get input and output MS name.
        // Those parameters were always called msin and msout.
        // However, SAS/MAC cannot handle a parameter and a group with the same
        // name, hence one can also use msin.name and msout.name.
        std::vector<string> inNames = parset.getStringVector ("msin.name",
                                                         std::vector<string>());
        if (inNames.empty()) {
          inNames = parset.getStringVector ("msin");
        }
        if (inNames.size() == 0)
          throw Exception("No input MeasurementSets given");
        // Find all file names matching a possibly wildcarded input name.
        // This is only possible if a single name is given.
        if (inNames.size() == 1) {
          if (inNames[0].find_first_of ("*?{['") != string::npos) {
            std::vector<string> names;
            names.reserve (80);
            casacore::Path path(inNames[0]);
            casacore::String dirName(path.dirName());
            casacore::Directory dir(dirName);
            // Use the basename as the file name pattern.
            casacore::DirectoryIterator dirIter (dir,
                                             casacore::Regex::fromPattern(path.baseName()));
            while (!dirIter.pastEnd()) {
              names.push_back (dirName + '/' + dirIter.name());
              dirIter++;
            }
            if (names.empty())
              throw Exception("No datasets found matching msin "
                       + inNames[0]);
            inNames = names;
          }
        }

        // Get the steps.
        // Currently the input MS must be given.
        // In the future it might be possible to have a simulation step instead.
        // Create MSReader step if input ms given.
        if (inNames.size() == 1) {
          reader = new MSReader (inNames[0], parset, "msin.");
        } else {
          reader = new MultiMSReader (inNames, parset, "msin.");
        }
        firstStep = DPStep::ShPtr (reader);
      }

      casacore::Path pathIn (reader->msName());
      casacore::String currentMSName (pathIn.absoluteName());

      // Create the other steps.
      std::vector<string> steps = parset.getStringVector (prefix + "steps");
      lastStep = firstStep;
      DPStep::ShPtr step;
      for (std::vector<string>::const_iterator iter = steps.begin();
           iter != steps.end(); ++iter) {
        string prefix(*iter + '.');
        // The alphabetic part of the name is the default step type.
        // This allows names like average1, out3.
        string defaulttype = (*iter);
        while (defaulttype.size()>0 && std::isdigit(*defaulttype.rbegin())) {
          defaulttype.resize(defaulttype.size()-1);
        }

        string type = parset.getString(prefix+"type", defaulttype);
        boost::algorithm::to_lower(type);
        // Define correct name for AOFlagger synonyms.
        if (type == "aoflagger" || type == "aoflag") {
          step = DPStep::ShPtr(new AOFlaggerStep (reader, parset, prefix));
        } else if (type == "averager"  ||  type == "average"  ||  type == "squash") {
          step = DPStep::ShPtr(new Averager (reader, parset, prefix));
        } else if (type == "madflagger"  ||  type == "madflag") {
          step = DPStep::ShPtr(new MedFlagger (reader, parset, prefix));
        } else if (type == "preflagger"  ||  type == "preflag") {
          step = DPStep::ShPtr(new PreFlagger (reader, parset, prefix));
        } else if (type == "uvwflagger"  ||  type == "uvwflag") {
          step = DPStep::ShPtr(new UVWFlagger (reader, parset, prefix));
        } else if (type == "counter"  ||  type == "count") {
          step = DPStep::ShPtr(new Counter (reader, parset, prefix));
        } else if (type == "phaseshifter"  ||  type == "phaseshift") {
          step = DPStep::ShPtr(new PhaseShift (reader, parset, prefix));
#ifdef HAVE_LOFAR_BEAM
        } else if (type == "demixer"  ||  type == "demix") {
          step = DPStep::ShPtr(new Demixer (reader, parset, prefix));
        } else if (type == "smartdemixer"  ||  type == "smartdemix") {
          step = DPStep::ShPtr(new DemixerNew (reader, parset, prefix));
        } else if (type == "applybeam") {
          step = DPStep::ShPtr(new ApplyBeam (reader, parset, prefix));
#endif
        } else if (type == "stationadder"  ||  type == "stationadd") {
          step = DPStep::ShPtr(new StationAdder (reader, parset, prefix));
        } else if (type == "scaledata") {
          step = DPStep::ShPtr(new ScaleData (reader, parset, prefix));
        } else if (type == "filter") {
          step = DPStep::ShPtr(new Filter (reader, parset, prefix));
        } else if (type == "applycal"  ||  type == "correct") {
          step = DPStep::ShPtr(new ApplyCal (reader, parset, prefix));
        } else if (type == "predict") {
          step = DPStep::ShPtr(new Predict (reader, parset, prefix));
        } else if (type == "h5parmpredict") {
          step = DPStep::ShPtr(new H5ParmPredict (reader, parset, prefix));
        } else if (type == "gaincal"  ||  type == "calibrate") {
          step = DPStep::ShPtr(new GainCal (reader, parset, prefix));
        } else if (type == "upsample") {
          step = DPStep::ShPtr(new Upsample (reader, parset, prefix));
        } else if (type == "split" || type == "explode") {
          step = DPStep::ShPtr(new Split (reader, parset, prefix));
        } else if (type == "ddecal") {
          step = DPStep::ShPtr(new DDECal (reader, parset, prefix));
        } else if (type == "interpolate") {
          step = DPStep::ShPtr(new Interpolate (reader, parset, prefix));
        } else if (type == "out" || type=="output" || type=="msout") {
          step = makeOutputStep(dynamic_cast<MSReader*>(reader), parset, prefix, currentMSName);
        } else {
          // Maybe the step is defined in a dynamic library.
          step = findStepCtor(type) (reader, parset, prefix);
        }
        if (lastStep) {
          lastStep->setNextStep (step);
        }
        lastStep = step;
        // Define as first step if not defined yet.
        if (!firstStep) {
          firstStep = step;
        }
      }
      // Add an output step if not explicitly added in steps (unless last step is a 'split' step)
      if (steps.size()==0 || (
          steps[steps.size()-1] != "out" &&
          steps[steps.size()-1] != "output" &&
          steps[steps.size()-1] != "msout" &&
          steps[steps.size()-1] != "split")) {
        step = makeOutputStep(dynamic_cast<MSReader*>(reader), parset, "msout.", currentMSName);
        lastStep->setNextStep (step);
        lastStep = step;
      }

      // Add a null step, so the last step can use getNextStep->process().
      DPStep::ShPtr nullStep(new NullStep());
      if (lastStep) {
        lastStep->setNextStep (nullStep);
      } else {
        firstStep = nullStep;
      }
      return firstStep;
    }

    DPStep::ShPtr DPRun::makeOutputStep (MSReader* reader,
                                         const ParameterSet& parset,
                                         const string& prefix,
                                         casacore::String& currentMSName)
    {
      DPStep::ShPtr step;
      casacore::String outName;
      bool doUpdate = false;
      if (prefix == "msout.") {
        // The last output step.
        outName = parset.getString ("msout.name", "");
        if (outName.empty()) {
          outName = parset.getString ("msout");
        }
      } else {
        // An intermediate output step.
        outName = parset.getString(prefix + "name");
      }

      // A name equal to . or the last name means an update of the last MS.
      if (outName.empty()  ||  outName == ".") {
        outName  = currentMSName;
        doUpdate = true;
      } else {
        casacore::Path pathOut(outName);
        if (currentMSName == pathOut.absoluteName()) {
          outName  = currentMSName;
          doUpdate = true;
        }
      }
      if (doUpdate) {
        // Create MSUpdater.
        // Take care the history is not written twice.
        // Note that if there is nothing to write, the updater won't do anything.
        step = DPStep::ShPtr(new MSUpdater(dynamic_cast<MSReader*>(reader),
                                           outName, parset, prefix,
                                           outName!=currentMSName));
      } else {
        step = DPStep::ShPtr(new MSWriter (reader, outName, parset, prefix));
        reader->setReadVisData (true);
      }
      casacore::Path pathOut(outName);
      currentMSName = pathOut.absoluteName();
      return step;
    }


  } //# end namespace
}
