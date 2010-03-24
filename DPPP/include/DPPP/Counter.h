//# Counter.h: DPPP step class to count flags
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

#ifndef DPPP_COUNTER_H
#define DPPP_COUNTER_H

// @file
// @brief DPPP step class to count flags

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/FlagCounter.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    // @ingroup NDPPP

    // This class is a DPStep class counting the number of flags per
    // baseline and channel.
    // It can be used for test purposes to know how many flags have been
    // set by the previous steps.

    class Counter: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Counter (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~Counter();

      // Process the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the average info.
      virtual void updateAverageInfo (AverageInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the flag counts.
      virtual void showCounts (std::ostream&) const;

    private:
      //# Data members.
      DPInput*    itsInput;
      string      itsName;
      uint        itsCount;
      FlagCounter itsFlagCounter;
    };

  } //# end namespace
}

#endif
