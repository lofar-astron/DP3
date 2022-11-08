// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Measurement Set constants which are not (yet) in casacore.
/// @author Lars Krombeen

#ifndef DPPP_MS_H
#define DPPP_MS_H

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

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

extern const std::string kAntennaTable;
extern const std::string kDataDescTable;
extern const std::string kObservationTable;

}  // namespace DP3MS

/**
 * Read the antenna set (e.g. LBA_OUTER, HBA_DUAL_INNER) from the OBSERVATION
 * table in the MS. This is a LOFAR-specific extension.
 * @param ms A measurement set, which may contain the antenna set field.
 * @return The antenna set of the MS, or an empty string if not found.
 */
std::string ReadAntennaSet(const casacore::MeasurementSet& ms);

}  // namespace base
}  // namespace dp3

#endif
