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

#ifndef __FLAGGER_FREQUENCYFLAGGER_H__
#define __FLAGGER_FREQUENCYFLAGGER_H__

#include <casa/Arrays.h>
#include <utility>
#include <vector>
#include <list>
#include <map>
#include <DPPP/Flagger.h>

/// @file
/// @brief Class to hold code for FrequencyFlagger in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    ///Foreward declarations
    class DataBuffer;
    class MsInfo;
    class RunDetails;
    class FlaggerStatistics;

    class FrequencyFlagger: public Flagger
    {
      public:
        FrequencyFlagger();
        ~FrequencyFlagger();

        /// All processing of one integration time happens in one go.
        void ProcessTimeslot(DataBuffer& data,
                             MsInfo& info,
                             RunDetails& details,
                             FlaggerStatistics& stats);

      protected:
      private:
        int FlagBaselineBand(casa::Matrix<casa::Bool>& Flags,
                             casa::Matrix<casa::Complex>& Data,
                             int flagCounter,
                             double FlagThreshold,
                             bool Existing,
                             int Algorithm);
        int NumChannels;
        int NumPolarizations;
    }; // FrequencyFlagger
  }; // namespace CS1
}; // namespace LOFAR

#endif //  __FLAGGER_FREQUENCYFLAGGER_H__
