//# FlagCounter.h: Class to keep counts of nr of flagged points
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

#ifndef DPPP_FLAGCOUNTER_H
#define DPPP_FLAGCOUNTER_H

// @file
// @brief Class to keep counts of nr of flagged points

#include <Common/lofar_vector.h>
#include <Common/lofar_string.h>
#include <Common/LofarTypes.h>
#include <casa/Arrays/Vector.h>

namespace LOFAR {
  namespace DPPP {

    // @ingroup NDPPP

    // This class contains counts the number of flags.
    // The flags can be counted per baseline, channel, and correlation.
    // Once the counting is completed, they can be printed using the 'show'
    // functions. When printing, the baselines counts are shown per antenna.

    class FlagCounter
    {
    public:
      // The default constructor creates an empty object.
      explicit FlagCounter (const string& comment)
        : itsComment (comment)
      {}

      // Setup 
      // Construct the object for the given phase direction and stations.
      void init (uint nbaselines, uint nchan, uint ncorr);

      // Increment the count per baseline.
      void incrBaseline (uint bl)
        { itsBLCounts[bl]++; }

      // Increment the count per channel.
      void incrChannel (uint chan)
        { itsChanCounts[chan]++; }

      // Increment the count per correlation.
      void incrCorrelation (uint corr)
        { itsCorrCounts[corr]++; }

      // Get the counts.
      const vector<int64>& baselineCounts() const
        { return itsBLCounts; }
      const vector<int64>& channelCounts() const
        { return itsChanCounts; }
      const vector<int64>& correlationCounts() const
        { return itsCorrCounts; }

      // Print the counts.
      void showBaseline (ostream& os, const casa::Vector<int>& ant1,
                         const casa::Vector<int>& ant2,
                         int64 npointsPerBaseline) const;
      void showChannel (ostream& os, int64 npointsPerChannel) const;
      void showCorrelation (ostream& os) const;

    private:
      vector<int64> itsBLCounts;
      vector<int64> itsChanCounts;
      vector<int64> itsCorrCounts;
      string        itsComment;
    };

  } //# end namespace
}

#endif
