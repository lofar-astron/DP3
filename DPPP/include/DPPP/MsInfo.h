//# Copyright (C) 2006-8
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
//# @author Adriaan Renting

#ifndef __CS1_PP_MS_INFO_H__
#define __CS1_PP_MS_INFO_H__

#include <ms/MeasurementSets.h>
#include <tables/Tables.h>
#include <utility>
#include <vector>
#include <string>
#include <map>

/// @file
/// @brief Class to hold code for retrieving information from MeasurementSet in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    class MsInfo
    {
    public:
      // Default constructor initializes to zero values.
      MsInfo();

      // Create from the MeasurementSet.
      // The orderedMainTable is the main table in time order. It can be the same
      // as the MeasurementSet if it is already in time order.
      MsInfo(const casa::MeasurementSet& ms, const casa::Table& orderedMainTable);

      ~MsInfo();

      // Get the baseline index of an antenna pair.
      // A negative value means that the baseline is not present.
      int getBaselineIndex (int ant1, int ant2) const
        { return BaselineIndex[ant1*NumAntennae + ant2]; }

      int                       NumSamples;
      int                       NumAntennae;
      int                       NumFields;
      int                       NumBands;
      int                       NumChannels;
      int                       NumPolarizations;
      int                       NumPairs;
      int                       NumTimeslots;
      double                    NoiseLevel;
      std::vector<casa::String> AntennaNames;
      casa::Vector<casa::Int>   Polarizations;
      double                    MaxBaselineLength;
      std::vector<int>          BaselineIndex;
      std::vector<double>       BaselineLengths;

      // Update the info for an output MS.
      void update (const casa::MeasurementSet& ms, int timestep);

      // print info to cout, for debugging.
      void                      PrintInfo(void);

    protected:
    private:
      // Calculate baseline length for baseline dependent flagging or filtering.
      void        ComputeBaselineLengths(const casa::MeasurementSet& MS);
    }; // class MsInfo
  }; // CS1
}; // namespace LOFAR

#endif // __CS1_PP_MS_INFO_H__
