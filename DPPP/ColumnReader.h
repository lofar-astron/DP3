// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step to read a column from the measurement set and overwrite the
/// buffer data with it.
/// @author Lars Krombeen

#ifndef DPPP_COLUMNREADER_H
#define DPPP_COLUMNREADER_H

#include "../Common/ParameterSet.h"
#include "DPStep.h"
#include "DPInput.h"

namespace DP3 {
namespace DPPP {

class ColumnReader : public DPStep {
 public:
  ColumnReader(DPInput& input, const ParameterSet&, const string& prefix,
               const string& column = "MODEL_DATA");

  /// Reads the given column for the measurementset for the rows of the input
  /// buffer, copies the buffer and change its data.
  bool process(const DPBuffer& buffer) override;

  /// Update the general info.
  void updateInfo(const DPInfo&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream& os, double duration) const override;

 private:
  DPInput& input_;           ///< Input MS to read the column from
  std::string name_;         ///< The name of the step (or prefix)
  std::string column_name_;  ///< Name of the column to use from the MS
  std::string operation_;    ///< Operation to use on the DATA column
  DPBuffer buffer_;          ///< Buffer to copy contents into
};

}  // namespace DPPP
}  // namespace DP3

#endif  // header guard
