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

#ifndef __FLAGGER_MADFLAGGER_H__
#define __FLAGGER_MADFLAGGER_H__

#include <casa/Arrays.h>
#include <utility>
#include <vector>
#include <list>
#include <map>
#include <DPPP/Flagger.h>

/// @file
/// @brief Class to hold code for MADFlagger in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    ///Forward declarations
    class DataBuffer;
    class MsInfo;
    class RunDetails;
    class FlaggerStatistics;

    // @ingroup IDPPP

    class MADFlagger: public Flagger
    {
      public:
        MADFlagger();
        ~MADFlagger();

        /// All processing of one integration time happens in one go.
        void ProcessTimeslot(DataBuffer& data,
                             MsInfo& info,
                             RunDetails& details,
                             FlaggerStatistics& stats);

      protected:
      private:
        casa::Cube<casa::Bool> MakeMask(const casa::Cube<casa::Float>& Values,
                                        const casa::Cube<casa::Bool>& Flags,
                                        double MaxLevel, int WindowSize,
                                        int Position, bool Exsisting);
        void ComputeThreshold(const casa::Cube<casa::Float>& Values,
                              const casa::Cube<casa::Bool>& Mask,
                              int TWindowSize, int FWindowSize,
                              int TimePos, int ChanPos, int PolPos,
                              float& Z1, float& Z2,
                              casa::Vector<casa::Float>& Medians);
        int FlagBaselineBand(casa::Cube<casa::Bool>& Flags,
                             const casa::Cube<casa::Float>& Data,
                             const casa::Cube<casa::Bool>& Mask,
                             int flagCounter,
                             double Threshold,
                             int Position, bool Existing,
                             int TWindowSize, int FWindowSize);
        int NumChannels;
        int NumPolarizations;
    }; // MADFlagger
  }; // CS1
}; // namespace LOFAR

#endif //  __FLAGGER_MADFLAGGER_H__
