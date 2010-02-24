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

#ifndef DPPP_UVWCALCULATOR_H
#define DPPP_UVWCALCULATOR_H

// @file
// @brief Class to calculate UVW coordinates

#include <Common/LofarLogger.h>
#include <Common/lofar_vector.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MBaseline.h>
#include <casa/Arrays/Vector.h>

namespace LOFAR {
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

      // Construct the object for the given phase direction and stations.
      UVWCalculator (const casa::MDirection& phaseDir,
                     const vector<casa::MPosition>& stationPositions);

      // Has the object been initialized?
      bool empty() const
        { return itsAntMB.empty(); }

      // get the UVW coordinates for the given baseline and time.
      casa::Vector<double> getUVW (uint ant1, uint ant2, double time);

    private:
      casa::MDirection              itsPhaseDir;
      casa::MeasFrame               itsFrame;
      vector<casa::MBaseline>       itsAntMB;
      vector<casa::Vector<double> > itsAntUvw;
      casa::Block<bool>             itsUvwFilled;
      double                        itsLastTime;
    };

  } //# end namespace
}

#endif
