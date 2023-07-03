// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_RESULTSTEP_H_
#define DP3_STEPS_RESULTSTEP_H_

#include <dp3/steps/Step.h>

namespace dp3 {
namespace steps {

/// @brief This class defines a step in the DP3 pipeline that keeps the result
/// to make it possible to get the result of another step.
/// It stores the result and does *NOT* call process() of the next step.

class ResultStep : public Step {
 public:
  /// Creates a MultiResultStep and sets a NullStep as its next step.
  ResultStep();

  ~ResultStep() override {}

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return {}; }

  /// Keep the buffer.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override {
    buffer_ = std::move(buffer);
    return true;
  }

  /// Finish does not do anything.
  void finish() override {}

  /// Show the step parameters.
  /// It does nothing.
  void show(std::ostream&) const override {}

  /// Get the result.
  /// Does not transfer ownership of the buffer to the caller. If that is
  /// required, use take() instead.
  const base::DPBuffer& get() const { return *buffer_; }

  /// Extract the result.
  /// Transfers ownership of the buffer to the caller.
  std::unique_ptr<base::DPBuffer> take() { return std::move(buffer_); }

 private:
  std::unique_ptr<base::DPBuffer> buffer_;
};

}  // namespace steps
}  // namespace dp3

#endif
