// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_COMBINE_H_
#define DP3_STEPS_COMBINE_H_

#include <dp3/steps/Step.h>

#include <xtensor/xtensor.hpp>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3::steps {

/// @brief DP3 step to combine two named buffers
/// This class is a DP3 step to combine to named buffers
class Combine final : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Combine(const common::ParameterSet&, const std::string& prefix);

  /// getRequiredFields should return all fields that process() reads.
  /// This implementation is merely an example.
  common::Fields getRequiredFields() const final { return kDataField; }

  /// getProvidedFields should return all fields that process() writes.
  /// This implementation is merely an example.
  common::Fields getProvidedFields() const final { return kDataField; }

  /// Process the data.
  bool process(std::unique_ptr<base::DPBuffer>) final;

  /// TODO: Process BDA data.
  bool process(std::unique_ptr<base::BdaBuffer>) final;

  /// Finish the processing of this step and subsequent steps.
  void finish() final;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) final;

  /// Set operation.
  void SetOperation(const std::string& operation);

  /// Show the step parameters.
  void show(std::ostream&) const final;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const final;

 private:
  enum class Operation { kAdd, kSubtract };

  std::string name_;
  common::NSTimer timer_;
  std::string buffer_name_;  //< Buffer that will be added or subtracted
  Operation operation_;      //< Add or subtract
};

}  // namespace dp3::steps

#endif
