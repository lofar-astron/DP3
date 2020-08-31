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
/// @author Lars Krombeen

#include "BDAMS.h"

namespace DP3 {
namespace DPPP {
namespace BDAMS {

/// BDA_TIME_AXIS table.
const std::string kBDATimeAxisTable = "BDA_TIME_AXIS";
const std::string kTimeAxisId = "BDA_TIME_AXIS_ID";
const std::string kIsBdaApplied = "IS_BDA_APPLIED";
const std::string kMaxTimeInterval = "MAX_TIME_INTERVAL";
const std::string kMinTimeInterval = "MIN_TIME_INTERVAL";
const std::string kUnitTimeInterval = "UNIT_TIME_INTERVAL";
const std::string kIntervalFactors = "INTEGER_INTERVAL_FACTORS";
const std::string kHasBDAOrdering = "HAS_BDA_ORDERING";
const std::string kFieldId = "FIELD_ID";
const std::string kSingleFactorPerBL = "SINGLE_FACTOR_PER_BASELINE";

/// BDA_FACTORS table.
const std::string kBDAFactorsTable = "BDA_FACTORS";
const std::string kFactor = "FACTOR";

/// SPECTRAL_WINDOW table.
const std::string kSpectralWindowTable = "SPECTRAL_WINDOW";
const std::string kSpectralWindowId = "SPECTRAL_WINDOW_ID";
const std::string kDataDescTable = "DATA_DESCRIPTION";
const std::string kBDAFreqAxisId = "BDA_FREQ_AXIS_ID";
const std::string kBDASetId = "BDA_SET_ID";
const std::string kChanFreq = "CHAN_FREQ";
const std::string kChanWidth = "CHAN_WIDTH";
const std::string kEffectiveBW = "EFFECTIVE_BW";
const std::string kResolution = "RESOLUTION";
const std::string kNumChan = "NUM_CHAN";
const std::string kTotalBandWidth = "TOTAL_BANDWIDTH";
const std::string kRefFrequency = "REF_FREQUENCY";

const std::string kLofarAntennaSet = "LOFAR_ANTENNA_SET";

const std::string kAntennaTable = "ANTENNA";
const std::string kName = "NAME";
const std::string kDishDiameter = "DISH_DIAMETER";
const std::string kPosition = "POSITION";
const std::string kMeasFreqRef = "MEAS_FREQ_REF";

const std::string kObservationTable = "OBSERVATION";

}  // namespace BDAMS
}  // namespace DPPP
}  // namespace DP3
