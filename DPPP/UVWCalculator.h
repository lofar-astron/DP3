//# UVWCalculator.h: Class to calculate UVW coordinates
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

// Note: this code is used by LOFAR and APERTIF software.

#ifndef DPPP_UVWCALCULATOR_H
#define DPPP_UVWCALCULATOR_H

// @file
// @brief Class to calculate UVW coordinates

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MBaseline.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MCBaseline.h>
#include <casacore/casa/Arrays/Vector.h>

namespace DP3 {
  namespace DPPP {

    // @ingroup NDPPP

    // This class calculates the UVW coordinates for a given baseline and
    // time stamp in the same way as done in LofarStMan.
    //
    // It calculates and caches the UVW coordinates per antenna and combines
    // them to get the baseline UVW coordinates. This is much faster than
    // calculating baseline UVW coordinates directly.

    class UVWCalculator
    {
    public:
      // The default constructor creates an empty object.
      UVWCalculator();

      // Construct the object for the given phase direction, array position,
      // and station positions.
      UVWCalculator (const casacore::MDirection& phaseDir,
                     const casacore::MPosition& arrayPosition,
                     const std::vector<casacore::MPosition>& stationPositions);

      // get the UVW coordinates for the given baseline and time.
      casacore::Vector<double> getUVW (uint ant1, uint ant2, double time);

    private:
      casacore::MDirection              itsPhaseDir;
      bool                          itsMovingPhaseDir;  
      casacore::MDirection::Convert     itsDirToJ2000;   //# direction to J2000
      casacore::MBaseline::Convert      itsBLToJ2000;    //# convert ITRF to J2000
      casacore::MeasFrame               itsFrame;
      std::vector<casacore::MBaseline>       itsAntMB;
      std::vector<casacore::Vector<double> > itsAntUvw;
      casacore::Block<bool>             itsUvwFilled;
      double                        itsLastTime;
    };

  } //# end namespace
}

#endif
