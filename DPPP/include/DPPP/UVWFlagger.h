//# UVWFlagger.h: DPPP step class to flag data on UVW coordinates
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

#ifndef DPPP_UVWFLAGGER_H
#define DPPP_UVWFLAGGER_H

// @file
// @brief DPPP step class to flag data on UVW coordinates

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/UVWCalculator.h>
#include <DPPP/FlagCounter.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  class ParameterSet;
  class ParameterValue;

  namespace DPPP {

    // @ingroup NDPPP

    // This class is a DPStep class flagging data points based on data
    // selections given in the parset file.
    // The following selections can be given:
    // <ul>
    //  <li> minimum and/or maximum UV distance
    //  <li> minimum or maximum value for U, V, and/or W
    //  <li> both can be used with a different phase center which can be
    //       be given as a position or as a moving source like SUN or JUPITER.
    // </ul>
    // The UVW values can be given in meters or in wavelengths.

    class UVWFlagger: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      // The antenna names are used to find antenna numbers.
      // The channel frequencies as they are in the input step must be given
      // starting at the start-channel.
      UVWFlagger (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~UVWFlagger();

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

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Test if uvw matches a range in meters.
      bool testUVWm (double uvw, const vector<double>& ranges);

      // Set flags for channels where uvw (in m) matches a range in wavelengths.
      void testUVWl (double uvw, const vector<double>& ranges,
                     bool* flagPtr, uint nrcorr);

      // Return a vector with UVW ranges.
      // It looks for the named parameter suffixed with 'range', 'min', and
      // 'max'. The returned vector contains 2 subsequent values for each range
      // (min and max are also turned into a range).
      // Optionally the values are squared to avoid having to take a sqrt
      // of the data's UVW coordinates.
      vector<double> fillUVW (const ParameterSet& parset,
                              const string& prefix,
                              const string& name,
                              bool square);

      // Handle the specification of a phase center.
      // It setups the UVWCalculator.
      void handleCenter();

      //# Data members.
      DPInput*             itsInput;
      string               itsName;
      uint                 itsNTimes;
      bool                 itsFlagOnWl; //# Flag on wavelengths?
      casa::Vector<double> itsRecWavel; //# reciprokes of wavelengths
      vector<double>       itsRangeUVm; //# UV ranges (in m) to be flagged
      vector<double>       itsRangeUm;  //# U  ranges (in m) to be flagged
      vector<double>       itsRangeVm;  //# V  ranges (in m) to be flagged
      vector<double>       itsRangeWm;  //# W  ranges (in m) to be flagged
      vector<double>       itsRangeUVl; //# UV ranges (in wavel) to be flagged
      vector<double>       itsRangeUl;  //# U  ranges (in wavel) to be flagged
      vector<double>       itsRangeVl;  //# V  ranges (in wavel) to be flagged
      vector<double>       itsRangeWl;  //# W  ranges (in wavel) to be flagged
      UVWCalculator        itsUVWCalc;
      vector<string>       itsCenter;
      NSTimer              itsTimer;
      NSTimer              itsUVWTimer;
    };

  } //# end namespace
}

#endif
