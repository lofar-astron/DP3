// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_ANTENNAFLAGGERSTEP_H
#define DP3_ANTENNAFLAGGERSTEP_H

#include <memory>
#include <string>
#include <vector>

#include <xtensor/xtensor.hpp>

#include "../antennaflagger/Flagger.h"

#include "InputStep.h"

namespace dp3 {
namespace steps {
class AntennaFlagger final : public Step {
 public:
  AntennaFlagger(InputStep* input, const common::ParameterSet& parset,
                 const string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField;
  }

  common::Fields getProvidedFields() const override { return kFlagsField; }

  void finish() override;
  void show(std::ostream& ostream) const override;
  bool process(const base::DPBuffer& buffer) override;
  void showTimings(std::ostream& ostream, double duration) const override;

 private:
  // Data statistics
  size_t n_antennas_;
  size_t n_channels_;
  size_t n_correlations_;

  // Data fields
  base::DPBuffer buffer_;
  xt::xtensor<bool, 2> selection_;

  std::string name_;
  std::string selection_string_;
  bool do_detect_outliers_;
  std::unique_ptr<dp3::antennaflagger::Flagger> flagger_;

  // AartfaacFlagger parameters
  float antenna_flagging_sigma_;
  size_t antenna_flagging_maxiters_;
  size_t antenna_flagging_outlier_threshold_;
  float station_flagging_sigma_;
  size_t station_flagging_maxiters_;

  // TImers
  common::NSTimer initialization_timer_;
  common::NSTimer computation_timer_;
  common::NSTimer flagging_timer_;
};

}  // namespace steps
}  // namespace dp3

#endif  // DP3_ANTENNAFLAGGERSTEP_H
