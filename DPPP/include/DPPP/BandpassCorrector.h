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

#ifndef __CORRECTOR_BANDPASSCORRECTOR_H__
#define __CORRECTOR_BANDPASSCORRECTOR_H__

#include <casa/Arrays.h>
#include <utility>
#include <vector>
#include <list>
#include <map>

/// @file
/// @brief Class to hold code for BandpassCorrector in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

/// This class applies a bandpass correction from a lookup table indexed with
/// RunDetails->Fixed

namespace LOFAR
{
  namespace CS1
  {
    //Forward declarations
    class MsInfo;
    class RunDetails;
    class DataBuffer;

    // @ingroup IDPPP

    class BandpassCorrector
    {
      public:
        BandpassCorrector(void);
        ~BandpassCorrector();
        /// All processing of one integration time happens in one go.
        void ProcessTimeslot(DataBuffer& data, MsInfo& info, RunDetails& details);
      protected:
      private:
        /// Processes one band of one baseline, fixed selects with profile is used
        void ProcessBaselineBand(casa::Matrix<casa::Complex>& In,
                                 casa::Matrix<casa::Complex>& Out,
                                 int fixed);
        int NumChannels; ///< From MsInfo
        int NumPolarizations; ///< From MsInfo
    }; // BandpassCorrector
  }; // namespace CS1
}; // namespace LOFAR

#endif //  __CORRECTOR_BANDPASSCORRECTOR_H__
