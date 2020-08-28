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

#ifndef DPPP_MS_H
#define DPPP_MS_H

#include <string>

using std::string;

namespace DP3 {
namespace DPPP {
namespace BDAMS {

/// BDA_TIME_AXIS table.
const string kBDATimeAxisTable = "BDA_TIME_AXIS";
const string kTimeAxisId = "BDA_TIME_AXIS_ID";
const string kIsBdaApplied = "IS_BDA_APPLIED";
const string kMaxTimeInterval = "MAX_TIME_INTERVAL";
const string kMinTimeInterval = "MIN_TIME_INTERVAL";
const string kUnitTimeInterval = "UNIT_TIME_INTERVAL";
const string kIntervalFactors = "INTEGER_INTERVAL_FACTORS";
const string kHasBDAOrdering = "HAS_BDA_ORDERING";
const string kFieldId = "FIELD_ID";
const string kSingleFactorPerBL = "SINGLE_FACTOR_PER_BASELINE";

/// BDA_FACTORS table.
const string kBDAFactorsTable = "BDA_FACTORS";
const string kFactor = "FACTOR";

/// SPECTRAL_WINDOW table.
const string kSpectralWindowTable = "SPECTRAL_WINDOW";
const string kSpectralWindowId = "SPECTRAL_WINDOW_ID";
const string kDataDescTable = "DATA_DESCRIPTION";
const string kBDAFreqAxisId = "BDA_FREQ_AXIS_ID";
const string kBDASetId = "BDA_SET_ID";
const string kChanFreq = "CHAN_FREQ";
const string kChanWidth = "CHAN_WIDTH";
const string kEffectiveBW = "EFFECTIVE_BW";
const string kResolution = "RESOLUTION";
const string kNumChan = "NUM_CHAN";
const string kTotalBandWidth = "TOTAL_BANDWIDTH";
const string kRefFrequency = "REF_FREQUENCY";

const string kLofarAntennaSet = "LOFAR_ANTENNA_SET";

const string kAntennaTable = "ANTENNA";
const string kName = "NAME";
const string kDishDiameter = "DISH_DIAMETER";
const string kPosition = "POSITION";
const string kMeasFreqRef = "MEAS_FREQ_REF";

}  // namespace BDAMS
}  // namespace DPPP
}  // namespace DP3

#endif