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

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <casa/Arrays/Cube.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    class Averager: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Averager (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~Averager();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the average info.
      virtual void updateAverageInfo (AverageInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

    private:
      // Average and return the result.
      DPBuffer average() const;

      // Copy the preAvg flags in the input buffer to the correct
      // time index in the output buffer.
      // If a flag is set, set all flags in corresponding PreAvg window.
      void copyPreAvgFlags (const casa::Cube<bool>& preAvgFlags,
                            const casa::Cube<bool>& flags,
                            int timeIndex);

      //# Data members.
      DPInput*        itsInput;
      string          itsName;
      DPBuffer        itsBuf;
      casa::Cube<int> itsNPoints;
      uint            itsNChanAvg;
      uint            itsNTimeAvg;
      uint            itsNTimes;
      double          itsTimeInterval;
    };

  } //# end namespace
}

#endif
