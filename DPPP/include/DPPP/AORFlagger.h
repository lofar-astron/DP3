//# AORFlagger.h: DPPP step class to flag data usibng rficonsole's functionality
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

#ifndef DPPP_AORFLAGGER_H
#define DPPP_AORFLAGGER_H

// @file
// @brief DPPP step class to flag using rficonsole's functionality

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/FlagCounter.h>
#include <Common/lofar_vector.h>
#include <Common/lofar_smartptr.h>
#include <AOFlagger/strategy/actions/strategyaction.h>
#include <AOFlagger/util/progresslistener.h>
#include <AOFlagger/quality/statisticscollection.h>

namespace LOFAR {

  namespace DPPP {
    class ParSet;

    // @ingroup NDPPP

    // This class is a DPStep class flagging data points based on the
    // rficonsole program (aka AOFlagger) written by Andre Offringa.
    // See the following papers for background information:
    // <ul>
    // <li> Post-correlation radio frequency interference classification
    //      methods -- http://arxiv.org/abs/1002.1957 
    // <li> A LOFAR RFI detection pipeline and its first results
    //      -- http://arxiv.org/abs/1007.2089
    //
    // When a correlation is flagged, all correlations for that data point
    // are flagged. It is possible to specify which correlations have to be
    // taken into account when flagging. Using, say, only XX may boost
    // performance with a factor 4, but miss points to be flagged.
    // It is also possible to specify the order in which the correlations
    // have to be tested.
    //
    // It is possible to flag specific baselines only using a selection on
    // baseline length.
    // <br>Furthermore it is possible to only flag the autocorrelations and
    // apply the resulting flags to the crosscorrelations, possibly selected
    // on baseline length.

    class AORFlagger: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      AORFlagger (DPInput*, const ParSet&, const string& prefix);

      virtual ~AORFlagger();

      // Process the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Write the statistics into the MS.
      virtual void addToMS (const string& msName);

      // Update the general info.
      // It is used to adjust the parms if needed.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the flagger counts.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Flag all baselines in the time window (using OpenMP to parallellize).
      // Process the buffers in the next step.
      void flag (uint rightOverlap);

      // Flag a single baseline using the rfistrategy.
      void flagBaseline (uint leftOverlap, uint windowSize,
                         uint rightOverlap, uint bl,
                         FlagCounter& counter,
			 rfiStrategy::Strategy&,
                         StatisticsCollection& rfiStats);

      // Add the flags to the statistics.
      void addStats (StatisticsCollection& rfiStats,
                     const Image2DPtr& reals, const Image2DPtr& imags,
		     const Mask2DCPtr& mask, const Mask2DPtr& origFlags,
                     int bl, uint polarization);

      // Fill the rfi strategy.
      void fillStrategy (boost::shared_ptr<rfiStrategy::Strategy>&);

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      uint             itsBufIndex;
      uint             itsNTimes;
      uint             itsNTimesToDo;
      string           itsStrategyName;
      uint             itsWindowSize;
      uint             itsOverlap;       //# extra time slots on both sides
      double           itsOverlapPerc;
      double           itsMemory;        //# Usable memory in GBytes
      double           itsMemoryPerc;
      double           itsMemoryNeeded;  //# Memory needed for data/flags
      bool             itsPulsarMode;
      bool             itsPedantic;
      bool             itsDoAutoCorr;
      bool             itsDoRfiStats;
      vector<DPBuffer> itsBuf;
      FlagCounter      itsFlagCounter;
      NSTimer          itsTimer;
      NSTimer          itsQualityTimer;  //# quality writing timer
      NSTimer          itsComputeTimer;  //# move/flag timer
      double           itsMoveTime;      //# data move timer (sum all threads)
      double           itsFlagTime;      //# flag timer (sum of all threads)
      double           itsQualTime;      //# quality timer (sum of all threads)
      boost::shared_ptr<rfiStrategy::Strategy> itsStrategy;
      DummyProgressListener itsProgressListener;
      StatisticsCollection  itsRfiStats;
      casa::Vector<double>  itsFreqs;
    };

  } //# end namespace
}

#endif
