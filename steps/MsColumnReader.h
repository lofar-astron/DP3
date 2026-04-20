// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step to read a column from the measurement set and overwrite the
/// buffer data with it.
/// @author Lars Krombeen

#ifndef DP3_STEPS_MSCOLUMNREADER_H_
#define DP3_STEPS_MSCOLUMNREADER_H_

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/Table.h>

#include <xtensor/containers/xtensor.hpp>

#include "steps/Step.h"

#include "common/ParameterSet.h"

namespace dp3 {
namespace steps {

class MsColumnReader final : public ModelDataStep {
 public:
  MsColumnReader(const common::ParameterSet&, const std::string& prefix,
                 MsType input_ms_type,
                 const std::string& column = "MODEL_DATA");

  common::Fields getRequiredFields() const override { return common::Fields(); }

  common::Fields getProvidedFields() const override { return kDataField; }

  /// Reads the given column for the measurementset for the rows of the input
  /// buffer and change its data.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  bool process(std::unique_ptr<base::BdaBuffer> buffer) override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream& os, double duration) const override;

  base::Direction GetFirstDirection() const override;

  bool accepts(MsType dt) const override { return dt == input_ms_type_; }

 private:
  casacore::Table table_;  ///< Input table to read the column from
  casacore::ArrayColumn<std::complex<float>> model_column_;
  casacore::Array<std::complex<float>> row_buffer_;
  std::string name_;         ///< The name of the step (or prefix)
  std::string column_name_;  ///< Name of the column to use from the MS
  MsType input_ms_type_;
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
