//# MedFlagger.h: DPPP step class to flag data based on median filtering
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

#ifndef DPPP_MEDFLAGGER_H
#define DPPP_MEDFLAGGER_H

// @file
// @brief DPPP step class to average in time and/or freq

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/FlagCounter.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    // @ingroup NDPPP

    // This class is a DPStep class flagging data points based on the median
    // of the absolute difference of the data and the median of the data.
    // Both medians are taken in a time/frequency window around the data point.
    // Only unflagged data points in the window are taken into account.
    // The size of the window is given in the parset file.
    // 
    // The window around data points at the edges is formed by mirroring the
    // data at the edge. For example, for channel 0 and a window size of 7
    // the data are mirrored, thus channels 3,2,1,0,1,2,3 will be used.
    // For channel 1 the channels 2,1,0,1,2,3,4 will be used.
    // The test program tMirror.cc can be used to check the correctness of
    // the alogorithm to determine the channels to use.
    //
    // Taking the median is an O(N) operation, thus doing it for all data
    // points is an O(N^2) operation. The test program tMedian.cc can be
    // used to test the performance of the algorithms to determine the median.
    // It shows that casacore's kthLargest outperforms STL's nth_element.
    // <br>
    // Shuffling the data around to be able to determine the medians is also
    // an expensive operation and takes as much time as the medians themselves.
    //
    // When a correlation is flagged, all correlations for that data point
    // are flagged. It is possible to specify which correlations have to be
    // taken into account when flagging. Using, say, only XX may boost
    // performance with a factor 4, but miss points to be flagged.
    // It is also possible to specify the order in which the correlations
    // have to be tested.

    class MedFlagger: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      MedFlagger (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~MedFlagger();

      // Process the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the average info.
      // It is used to adjust the parms if needed.
      virtual void updateAverageInfo (AverageInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the flagger counts.
      virtual void showCounts (std::ostream&) const;

      // Flag for the entry at the given index.
      // Use the given time entries for the medians.
      // Process the result in the next step.
      void flag (uint index, const vector<uint>& timeEntries);

      // Compute the median factors for given baseline, channel, and
      // correlation.
      void computeFactors (const vector<uint>& timeEntries,
                           uint bl, int chan, int corr,
                           int nchan, int ncorr,
                           float& Z1, float& Z2,
                           float* tempBuf);

    private:
      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      float            itsThreshold;
      uint             itsFreqWindow;
      uint             itsTimeWindow;
      uint             itsNTimes;
      uint             itsNTimesDone;
      vector<uint>     itsFlagCorr;
      vector<DPBuffer> itsBuf;
      FlagCounter      itsFlagCounter;
    };

  } //# end namespace
}

#endif
