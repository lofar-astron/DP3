// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step to read a column from the measurement set and overwrite the
/// buffer data with it.
/// @author Lars Krombeen

#ifndef DP3_STEPS_COLUMNREADER_H_
#define DP3_STEPS_COLUMNREADER_H_

#include "../common/ParameterSet.h"
#include <dp3/steps/Step.h>
#include "InputStep.h"

namespace dp3 {
namespace steps {

class ColumnReader : public ModelDataStep {
 public:
  ColumnReader(InputStep& input, const common::ParameterSet&,
               const string& prefix, const string& column = "MODEL_DATA");

  common::Fields getRequiredFields() const override {
    common::Fields fields;
    if ((operation_ == Operation::kAdd) ||
        (operation_ == Operation::kSubtract)) {
      fields |= kDataField;
    }
    return fields;
  }

  common::Fields getProvidedFields() const override { return kDataField; }

  /// Reads the given column for the measurementset for the rows of the input
  /// buffer, copies the buffer and change its data.
  bool process(const base::DPBuffer& buffer) override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream& os, double duration) const override;

  base::Direction GetFirstDirection() const override;

 private:
  enum class Operation { kReplace, kAdd, kSubtract };

  InputStep& input_;         ///< Input MS to read the column from
  std::string name_;         ///< The name of the step (or prefix)
  std::string column_name_;  ///< Name of the column to use from the MS
  Operation operation_;      ///< Operation to use on the DATA column
  base::DPBuffer buffer_;    ///< Buffer to copy contents into
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
