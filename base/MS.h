// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Measurement Set constants which are not (yet) in casacore.
/// @author Lars Krombeen

#ifndef DPPP_MS_H
#define DPPP_MS_H

#include <string>

namespace dp3 {
namespace base {
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
extern const std::string kSpectralWindowId;

/// SPECTRAL_WINDOW table.
extern const std::string kSpectralWindowTable;
extern const std::string kBDAFreqAxisId;
extern const std::string kBDASetId;

extern const std::string kLofarAntennaSet;

extern const std::string kAntennaTable;
extern const std::string kDataDescTable;
extern const std::string kObservationTable;

}  // namespace DP3MS
}  // namespace base
}  // namespace dp3

#endif
