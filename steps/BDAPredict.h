// BDAPredict.h: class to directly predict baseline dependent averaged (BDA)
// visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief class to directly predict baseline dependent averaged (BDA)
/// visibilities from a source model
/// @author Sebastiaan van der Tol

#ifndef DP3_BDAPREDICT_H
#define DP3_BDAPREDICT_H

#include "InputStep.h"

#include "../base/BDABuffer.h"

#include <queue>
#include <utility>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DP3 step class to predict BDA visibilities from a source model
/// @author Sebastiaan van der Tol

class BDAPredict : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  BDAPredict(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~BDAPredict();

  /// Process the data.
  /// Incoming BDABuffers are buffered in a queue
  /// and send to the the next step when all baseline groups are complete.
  /// This is necessary because baseline groups may overlap multiple BDABuffers,
  /// while the predict is done by calls to the (regular) Predict step, which
  /// needs the baseline groups to be complete
  bool process(std::unique_ptr<base::BDABuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Boolean if this step can process this type of data.
  bool accepts(MSType dt) const override { return dt == BDA; }

  /// Return which datatype this step outputs.
  MSType outputs() const override { return BDA; }

 private:
  InputStep* input_;

  // Need to store a reference to the parset to create the Predict substeps in
  // updateInfo()
  const common::ParameterSet& parset_;
  std::string name_;

  // Structure to keep track how many rows in a BDABuffer have been processed
  struct BufferInfo {
    std::unique_ptr<base::BDABuffer> buffer;
    size_t nr_rows_filled;
  };

  /// queue buffering the incoming BDABuffers
  std::queue<BufferInfo> buffers_;

  /// class representing a group of baselines with the same averaging parameters
  class BaselineGroup;

  /// map the pair of averaging parameters (time and frequency) to baseline
  /// group
  std::map<std::pair<int, int>, BaselineGroup> averaging_to_baseline_group_map_;

  /// vector, interpreted as a map of baseline index to
  /// a pair consisting of the baseline group it belongs to and index within
  /// that group
  std::vector<std::pair<BaselineGroup*, int>> index_to_baseline_group_map_;

  common::NSTimer timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
