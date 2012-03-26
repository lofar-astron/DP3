//# StationAdder.h: DPPP step class to add station to a superstation
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
// @brief DPPP step class to average in time and/or freq

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <Common/ParameterRecord.h>

namespace LOFAR {

  namespace DPPP {
    class ParSet;

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
      StationAdder (DPInput*, const ParSet&, const string& prefix);

      virtual ~StationAdder();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      //# Data members.
      DPInput*        itsInput;
      string          itsName;
      DPBuffer        itsBuf;
      ParameterRecord itsStatRec;     // stations definitions
      vector<int>     itsStations;    // >=0: superstation id of each station
      vector<string>  itsNewNames;    // Names of new superstation
      uint            itsMinNStation; // flag data if too few unflagged stations
      bool            itsUseWeight;   // False = use weight 1 per station
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
