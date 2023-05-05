// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_DUMMYSTEP_H_
#define DP3_STEPS_DUMMYSTEP_H_

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step class that does nothing.
/// This class is an empty Step subclass to use as implementation template.
class DummyStep : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DummyStep(const common::ParameterSet&, const std::string& prefix);

  /// getRequiredFields should return all fields that process() reads.
  /// This implementation is merely an example.
  common::Fields getRequiredFields() const override {
    return kWeightsField | kUvwField;
  }

  /// getProvidedFields should return all fields that process() writes.
  /// This implementation is merely an example.
  common::Fields getProvidedFields() const override { return kWeightsField; }

  /// Process the data. The dummy step forwards the data to its next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  /// Process BDA data. The dummy step forwards the data to its next step.
  bool process(std::unique_ptr<base::BDABuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

 private:
  std::string name_;
  common::NSTimer timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
