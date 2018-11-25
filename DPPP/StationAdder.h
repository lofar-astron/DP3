//# StationAdder.h: DPPP step class to add stations as a superstation
//# Copyright (C) 2012
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

#ifndef DPPP_STATIONADDER_H
#define DPPP_STATIONADDER_H

// @file
// @brief DPPP step class to add stations as a superstation

#include "DPInput.h"
#include "DPBuffer.h"
#include "UVWCalculator.h"

#include "../Common/ParameterRecord.h"

#include <casacore/measures/Measures/MPosition.h>

namespace DP3 {
  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class summing stations to a superstation.
    //
    // It is possible to define one or more groups of stations to be summed.
    // Each group has a name which is the name of the new station.
    // The complex values of baselines are added for which one station occurs
    // in only one group. A baseline is not added if no or both stations are
    // member of a summing group.
    // <br>The summation is done in a weighted way, where the weight of a
    // new station is the sum of the original weights. Optionally weights 1
    // can be used instead of the original weights.
    //
    // Only unflagged data points are used. If too few data points are
    // unflagged, the output data point is flagged.
    //
    // Questions:
    // 1. check if phases do not differ too much? Flag if too much?
    // 2. must all stations exist or possible that some don't?

    class StationAdder: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      StationAdder (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~StationAdder();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Add new meta info to the MS.
      virtual void addToMS (const string& msName);

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      // Return the indices of the stations in antennaNames matching
      // the pattern list.
      // The patterns are processed from left to right. A pattern can start
      // with ! or ^ meaning that the the matches are discarded. In this
      // way first a broad pattern can be given, which can be narrowed down.
      // A warning is given if a pattern does not match any station name.
      static std::vector<int> getMatchingStations
      (const casacore::Vector<casacore::String>& antennaNames,
       const std::vector<string>& patterns);

    private:
      // Update the beam info subtables.
      void updateBeamInfo (const string& msName, uint origNant,
                           casacore::Table& antTab);

      //# Data members.
      DPInput*        itsInput;
      string          itsName;
      DPBuffer        itsBuf;
      DPBuffer        itsBufTmp;
      ParameterRecord itsStatRec;     // stations definitions
      std::vector<casacore::Vector<int> > itsParts;  // the stations in each superstation
      std::vector<std::vector<int> > itsBufRows; // old baseline rows in each new baseline
      uint            itsMinNPoint  ;  // flag data if too few unflagged data
      bool            itsMakeAutoCorr; // also form new auto-correlations?
      bool            itsSumAutoCorr;  // sum auto- or cross-correlations?
      bool            itsDoAverage;    // average or sum?
      bool            itsUseWeight;    // false = use weight 1 per station
      UVWCalculator   itsUVWCalc;
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
