// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_COMBINE_H_
#define DP3_STEPS_COMBINE_H_

#include "Step.h"

#include <xtensor/containers/xtensor.hpp>

#include "common/ParameterSet.h"
#include "common/Timer.h"

namespace dp3::steps {

/// @brief DP3 step to combine two named buffers
/// This class is a DP3 step to combine to named buffers
class Combine final : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// @p ms_type specifies if BDA or regular data is combined.
  Combine(const common::ParameterSet&, const std::string& prefix,
          MsType ms_type);

  /// getRequiredFields should return all fields that process() reads.
  /// This implementation is merely an example.
  common::Fields getRequiredFields() const final { return kDataField; }

  /// getProvidedFields should return all fields that process() writes.
  /// This implementation is merely an example.
  common::Fields getProvidedFields() const final { return kDataField; }

  bool process(std::unique_ptr<base::DPBuffer>) final;

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

  bool accepts(MsType data_type) const override {
    return data_type == ms_type_;
  }

  MsType outputs() const override { return ms_type_; }

 private:
  enum class Operation { kAdd, kSubtract };

  std::string name_;
  common::NSTimer timer_;
  std::string buffer_name_;  //< Buffer that will be added or subtracted
  Operation operation_;      //< Add or subtract
  MsType ms_type_;           //< BDA or regular data
};

}  // namespace dp3::steps

#endif
