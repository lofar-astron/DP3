//# MedFlagger.h: DPPP step class to average in time and/or freq
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

#include <DPPP/DPStep.h>
#include <DPPP/DPBuffer.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    class MedFlagger: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      explicit MedFlagger (const ParameterSet&, const string& prefix);

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
      virtual void show (std::ostream&);

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
      string           itsName;
      float            itsThreshold;
      uint             itsFreqWindow;
      uint             itsTimeWindow;
      uint             itsNTimes;
      uint             itsNTimesDone;
      vector<uint>     itsFlagCorr;
      vector<DPBuffer> itsBuf;
    };

  } //# end namespace
}

#endif
