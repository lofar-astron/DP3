// BdaGroupPredict.h: class to directly predict baseline dependent averaged
// (BDA) visibilities from a source model Copyright (C) 2020 ASTRON (Netherlands
// Institute for Radio Astronomy) SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief class to directly predict baseline dependent averaged (BDA)
/// visibilities from a source model
/// @author Sebastiaan van der Tol

#ifndef DP3_BDAGROUPPREDICT_H
#define DP3_BDAGROUPPREDICT_H

#include "InputStep.h"

#include "base/BdaBuffer.h"

#include <map>
#include <queue>
#include <utility>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DP3 step class to predict BDA visibilities from a source model
/// @author Sebastiaan van der Tol

class BdaGroupPredict : public ModelDataStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  BdaGroupPredict(const common::ParameterSet&, const std::string& prefix);

  /**
   * Constructs the object with explicit source patterns.
   * @param input_step Input step, for reading extra data.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param source_patterns Source patterns.
   */
  BdaGroupPredict(const common::ParameterSet&, const std::string& prefix,
                  const std::vector<std::string>& source_patterns);

  ~BdaGroupPredict() override;

  common::Fields getRequiredFields() const override;

  /// Processes the data.
  /// Buffers incoming BdaBuffers in a queue and sends them to the the next step
  /// when all baseline groups are complete.
  /// This is necessary because baseline groups may overlap multiple BdaBuffers,
  /// while the predict is done by calls to the Predict step, which needs
  /// complete baseline groups.
  bool process(std::unique_ptr<base::BdaBuffer>) override;

  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

  MsType outputs() const override { return MsType::kBda; }

  base::Direction GetFirstDirection() const override;

 private:
  // Need to store a reference to the parset to create the OnePredict
  // substeps in updateInfo()
  const common::ParameterSet& parset_;
  std::string name_;

  // Structure to keep track how many rows in a BdaBuffer have been processed
  struct BufferInfo {
    std::unique_ptr<base::BdaBuffer> buffer;
    size_t nr_rows_filled;
  };

  /// queue buffering the incoming BdaBuffers
  std::queue<BufferInfo> buffers_;

  std::vector<std::string> source_patterns_;

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
