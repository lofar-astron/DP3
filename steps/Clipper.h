// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_CLIPPER_H_
#define DP3_STEPS_CLIPPER_H_

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "OnePredict.h"
#include "PreFlagger.h"
#include "ResultStep.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step for clipping bright sources from the data
/// This class is a Step class that predicts visibilities from
/// a skymodel and clips the bright sources from the input data.

class Clipper : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Clipper(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kFlagsField | kUvwField;
  }

  common::Fields getProvidedFields() const override { return kFlagsField; }

  /// Process the data. The step forwards the data to its next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

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
  size_t counter_;
  size_t time_step_;
  float ampl_max_;
  std::shared_ptr<OnePredict> predict_step_;
  std::shared_ptr<ResultStep> result_step_;
  std::unique_ptr<dp3::base::DPBuffer> predict_buffer_;
  base::DPBuffer::FlagsType last_flags_;
};

}  // namespace steps
}  // namespace dp3

#endif
