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

    //# Forward Declarations.
    class DPInfo;
    class ParSet;

    // @ingroup NDPPP

    // This class contains counts the number of flags.
    // The flags can be counted per baseline, channel, and correlation.
    // Once the counting is completed, they can be printed using the 'show'
    // functions. When printing, the baselines counts are shown per antenna.
    //
    // Optionally the flagging percentages can be saved in a table.
    // The name of the table is the MS name suffixed by the step name and '.flagxx'.

    class FlagCounter
    {
    public:
      // The default constructor creates an emty object. It does not save.
      FlagCounter();

      // The constructor creates an empty object.
      // It reads info from the parset to see if percentages have to be saved.
      FlagCounter (const string& msName, const ParSet&, const string& prefix);

      // Size all counters and initialize them to zero using the sizes
      // from the DPInfo object.
      void init (const DPInfo& info);

      // Size all counters to that's sizes and initialize them to zero.
      ///      void init (const FlagCounter& that);

      // Increment the count per baseline.
      void incrBaseline (uint bl)
        { itsBLCounts[bl]++; }

      // Increment the count per channel.
      void incrChannel (uint chan)
        { itsChanCounts[chan]++; }

      // Increment the count per correlation.
      void incrCorrelation (uint corr)
        { itsCorrCounts[corr]++; }

      // Add the contents of that to this.
      void add (const FlagCounter& that);

      // Get the counts.
      const vector<int64>& baselineCounts() const
        { return itsBLCounts; }
      const vector<int64>& channelCounts() const
        { return itsChanCounts; }
      const vector<int64>& correlationCounts() const
        { return itsCorrCounts; }

      // Print the counts and optionally save percentages in a table.
      void showBaseline    (ostream& os, int64 ntimes) const;
      void showChannel     (ostream& os, int64 ntimes) const;
      void showCorrelation (ostream& os, int64 ntimes) const;

      // Show percentage with 1 decimal.
      static void showPerc1 (std::ostream&, double value, double total);

      // Show percentage with 3 decimals.
      static void showPerc3 (std::ostream&, double value, double total);

    private:
      // Save the percentages per station in a table.
      void saveStation (int64 npoints, const casa::Vector<int64>& nused,
                        const casa::Vector<int64>& count) const;

      // Save the percentages per channel.
      void saveChannel (int64 npoints,
                        const casa::Vector<int64>& count) const;

      //# Data members.
      const DPInfo* itsInfo;
      string        itsSaveName;
      double        itsWarnPerc;
      bool          itsShowFF;
      vector<int64> itsBLCounts;
      vector<int64> itsChanCounts;
      vector<int64> itsCorrCounts;
    };

  } //# end namespace
}

#endif
