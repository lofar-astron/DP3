// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_ANTENNAFLAGGER_H_
#define DP3_STEPS_ANTENNAFLAGGER_H_

#include <memory>
#include <string>

#include "Step.h"

#include "../antennaflagger/Flagger.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {
class AntennaFlagger final : public Step {
 public:
  explicit AntennaFlagger(const common::ParameterSet& parset,
                          const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField;
  }

  common::Fields getProvidedFields() const override { return kFlagsField; }

  void updateInfo(const base::DPInfo& info) override;
  void finish() override;
  void show(std::ostream& ostream) const override;
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;
  void showTimings(std::ostream& ostream, double duration) const override;

 private:
  std::string name_;
  common::BaselineOrder baseline_order_;
  std::unique_ptr<dp3::antennaflagger::Flagger> flagger_;

  // AartfaacFlagger parameters
  float antenna_flagging_sigma_;
  size_t antenna_flagging_max_iterations_;
  size_t antenna_flagging_outlier_threshold_;
  float station_flagging_sigma_;
  size_t station_flagging_max_iterations_;

  // TImers
  common::NSTimer initialization_timer_;
  common::NSTimer computation_timer_;
  common::NSTimer flagging_timer_;
};

}  // namespace steps
}  // namespace dp3

#endif  // DP3_ANTENNAFLAGGERSTEP_H
