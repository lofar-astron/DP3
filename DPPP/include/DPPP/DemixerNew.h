//# DemixerNew.h: DPPP step class to subtract A-team sources in adaptive way
//# Copyright (C) 2013
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
//# $Id: DemixerNew.h 23223 2012-12-07 14:09:42Z schoenmakers $
//#
//# @author Ger van Diepen

#ifndef DPPP_DEMIXERNEW_H
#define DPPP_DEMIXERNEW_H

// @file
// @brief DPPP step class to subtract A-team sources in adaptive way

#include <DPPP/DemixInfo.h>
#include <DPPP/DemixWorker.h>
#include <DPPP/DPInput.h>
#include <DPPP/Filter.h>
#include <ParmDB/ParmDB.h>
#include <Common/lofar_smartptr.h>
#include <Common/lofar_map.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to subtract the strong A-team sources
    // in a smart way (algorithm developed by Reinout van Weeren).
    // It operates as follows:
    // <ol>
    // <li> Per time window (default 2 min) demixing is done separately.
    // <li> Using a simple Ateam model (the A and other strong sources)
    //      and the beam model, the StokesI flux is predicted per source.
    //      The antennae of baselines with flux>threshold are counted.
    //      Only an antenna counted in at least min_antenna baselines is
    //      solved for. If not such antennae exist, the source is ignored.
    // <li> The target is predicted as well. Note that an Ateam source within
    //      the target field is removed from the Ateam model and added/replaced
    //      in the target model.
    //      The target model is usually created using gsm.py.
    // <li> For the core baselines the ratio target/Ateam flux is used to
    //      determine if the target has to be included, ignored or deprojected.
    // </ol>
    // It is based on the demixing.py script made by Bas vd Tol and operates
    // per time chunk as follows:
    // <ul>
    //  <li> The data are phase-shifted and averaged for each source.
    //  <li> Demixing is done using the combined results.
    //  <li> For each source a BBS solve, smooth, and predict is done.
    //  <li> The predicted results are subtracted from the averaged data.
    // </ul>

    class DemixerNew: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      DemixerNew (DPInput*, const ParameterSet&, const string& prefix);

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the counts.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Process the data collected in itsBuf.
      void processData();

      // Export the solutions to a ParmDB.
      void writeSolutions (double startTime, int ntime);

      // Add the mean and M2 (square of differences) of a part in a
      // numerically stable way.
      void addMeanM2 (double& mean, double& m2, size_t& nr,
                      double partmean, double partm2, size_t partnr) const;

      // Show a statistic.
      void showStat (ostream& os, double n, double ntot,
                     const string& str1, const string& str2) const;

      // Show a percentage with 1 decimal.
      void showPerc1 (ostream& os, float perc) const;

      //# Data members.
      DPInput*                itsInput;
      string                  itsName;
      DemixInfo               itsDemixInfo;
      string                  itsInstrumentName;
      shared_ptr<BBS::ParmDB> itsParmDB;
      Filter                  itsFilter;    //# only used for getInfo()
      vector<DemixWorker>     itsWorkers;
      vector<DPBuffer>        itsBufIn;
      vector<DPBuffer>        itsBufOut;
      vector<vector<double> > itsSolutions; //# all solutions in a time window
      map<string,int>         itsParmIdMap; //# -1 = new parm name
      uint                    itsNTime;
      uint                    itsNTimeOut;
      uint                    itsNChunk;
      //# Timers.
      NSTimer itsTimer;
      NSTimer itsTimerDemix;
      NSTimer itsTimerDump;  //# writeSolutions
      NSTimer itsTimerNext;  //# next step (parallel to writeSolutions)
    };

  } //# end namespace
} //# end namespace

#endif
