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

 private:
  DPInput& input_;           ///< Input MS to read the column from
  std::string name_;         ///< The name of the step (or prefix)
  std::string column_name_;  ///< Name of the column to use from the MS
};

}  // namespace DPPP
}  // namespace DP3

#endif  // header guard
