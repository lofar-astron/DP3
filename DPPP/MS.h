// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Measurement Set constants which are not (yet) in casacore.
/// @author Lars Krombeen

#ifndef DPPP_MS_H
#define DPPP_MS_H

#include <string>

namespace DP3 {
namespace DPPP {
namespace DP3MS {  // Avoid name conflict with casacore::MS.

/// BDA_TIME_AXIS table.
extern const std::string kBDATimeAxisTable;
extern const std::string kTimeAxisId;
extern const std::string kIsBdaApplied;
extern const std::string kMaxTimeInterval;
extern const std::string kMinTimeInterval;
extern const std::string kUnitTimeInterval;
extern const std::string kIntervalFactors;
extern const std::string kHasBDAOrdering;
extern const std::string kFieldId;
extern const std::string kSingleFactorPerBL;

/// BDA_FACTORS table.
extern const std::string kBDAFactorsTable;
extern const std::string kFactor;

/// SPECTRAL_WINDOW table.
extern const std::string kSpectralWindowTable;
extern const std::string kSpectralWindowId;
extern const std::string kDataDescTable;
extern const std::string kBDAFreqAxisId;
extern const std::string kBDASetId;
extern const std::string kChanFreq;
extern const std::string kChanWidth;
extern const std::string kEffectiveBW;
extern const std::string kResolution;
extern const std::string kNumChan;
extern const std::string kTotalBandWidth;
extern const std::string kRefFrequency;

extern const std::string kLofarAntennaSet;

extern const std::string kAntennaTable;
extern const std::string kName;
extern const std::string kDishDiameter;
extern const std::string kPosition;
extern const std::string kMeasFreqRef;

extern const std::string kObservationTable;

}  // namespace DP3MS
}  // namespace DPPP
}  // namespace DP3

#endif
