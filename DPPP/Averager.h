//# Averager.h: DPPP step class to average in time and/or freq
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

#ifndef DPPP_AVERAGER_H
#define DPPP_AVERAGER_H

// @file
// @brief DPPP step class to average in time and/or freq

#include "DPInput.h"
#include "DPBuffer.h"

#include <casacore/casa/Arrays/Cube.h>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class calculating the weighted average of
    // data in time and/or frequency.
    // <br>
    // Only unflagged data points are used. The average is calculated as
    // <tt>sum(data*weight) / sum(weight)</tt> and the sum of the weights
    // is the weight of the new data point. If all data point to use are
    // flagged, the resulting data point and weight are set to zero and flagged.
    //
    // It keeps track of the FullResFlags. It sets them if the corresponding
    // data point is flagged. Note that multiple FullResFlags elements map to
    // a single data point if some averaging was done before.

    class Averager: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Averager (DPInput*, const ParameterSet&, const string& prefix);

      // Construct the object using the given parameters.
      Averager (DPInput*, const string& stepname,
                uint nchanAvg, uint ntimeAvg);

      virtual ~Averager();

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

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Average into itsBufOut.
      void average();

      // Copy the fullRes flags in the input buffer to the correct
      // time index in the output buffer.
      // If a flag is set, set all flags in corresponding FullRes window.
      void copyFullResFlags (const casacore::Cube<bool>& fullResFlags,
                             const casacore::Cube<bool>& flags,
                             int timeIndex);

      // Get the value in Hertz of a string like "1000 MHz". If unit is
      // omitted it defaults to Hertz
      double getFreqHz(const string& freqstr);

      //# Data members.
      DPInput*        itsInput;
      string          itsName;
      DPBuffer        itsBuf;
      DPBuffer        itsBufTmp;
      DPBuffer        itsBufOut;
      casacore::Cube<int> itsNPoints;
      casacore::Cube<casacore::Complex> itsAvgAll;
      casacore::Cube<float>         itsWeightAll;
      casacore::Cube<bool>          itsFullResFlags;
      double          itsFreqResolution;
      double          itsTimeResolution;
      uint            itsNChanAvg;
      uint            itsNTimeAvg;
      uint            itsMinNPoint;
      float           itsMinPerc;
      uint            itsNTimes;
      double          itsTimeInterval;
      bool            itsNoAvg;           //# No averaging (i.e. both 1)?
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
