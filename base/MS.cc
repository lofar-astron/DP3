// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @author Lars Krombeen

#include "MS.h"

#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableRecord.h>

namespace dp3 {
namespace base {
namespace DP3MS {

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
const std::string kSpectralWindowId = "SPECTRAL_WINDOW_ID";

/// SPECTRAL_WINDOW table extensions.
const std::string kSpectralWindowTable = "SPECTRAL_WINDOW";
const std::string kBDAFreqAxisId = "BDA_FREQ_AXIS_ID";
const std::string kBDASetId = "BDA_SET_ID";

const std::string kAntennaTable = "ANTENNA";
const std::string kDataDescTable = "DATA_DESCRIPTION";
const std::string kObservationTable = "OBSERVATION";

}  // namespace DP3MS

std::string ReadAntennaSet(const casacore::MeasurementSet& ms) {
  const std::string kLofarAntennaSet = "LOFAR_ANTENNA_SET";
  const casacore::Table observation_table(
      ms.keywordSet().asTable(DP3MS::kObservationTable));
  std::string antenna_set;
  if (observation_table.nrow() > 0 &&
      observation_table.tableDesc().isColumn(kLofarAntennaSet)) {
    antenna_set = casacore::ScalarColumn<casacore::String>(observation_table,
                                                           kLofarAntennaSet)(0);
  }
  return antenna_set;
}

}  // namespace base
}  // namespace dp3
