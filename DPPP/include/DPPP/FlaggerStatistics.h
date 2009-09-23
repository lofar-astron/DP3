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

#ifndef __CS1_PP_FLAGGER_STATISTICS_H__
#define __CS1_PP_FLAGGER_STATISTICS_H__

#include <casa/Arrays.h>
#include <iostream>
#include <vector>
#include <DPPP/MsInfo.h>

/// @file
/// @brief Class to hold code for holding the statistics of a Flagger pass in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    class FlaggerStatistics
    {
      public:
        FlaggerStatistics(MsInfo& info);
        ~FlaggerStatistics();
        /// Will output formatted statistics to the output stream (usually cout)
        void PrintStatistics(std::ostream& output);
        int& operator()(int x, int y, int z); ///< for quick indexing of the internal data

      protected:
      private:
        int                       NumAntennae;
        int                       NumBands;
        casa::Cube< int >         Statistics; ///< A cube of antenna x antenna x bands
        int                       Normalizer; ///< the total count of antenna x antenna x bands
        std::vector<casa::String> AntennaNames;
    }; // FlaggerStatistics
  }; // CS1
}; // namespace LOFAR

#endif //  __CS1_PP_FLAGGER_STATISTICS_H__
