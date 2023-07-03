// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_MULTIRESULTSTEP_H_
#define DP3_STEPS_MULTIRESULTSTEP_H_

#include <dp3/steps/Step.h>

namespace dp3 {
namespace steps {

/// @brief This class defines step in the DP3 pipeline that keeps the result
/// to make it possible to get the result of another step.
/// It keeps the result and calls process of the next step.
/// Buffers are accumulated until cleared.
class MultiResultStep : public Step {
 public:
  /// Creates a MultiResultStep and sets a NullStep as its next step.
  /// @param size The number of buffers the MultiResultStep should store.
  explicit MultiResultStep(unsigned int size);

  ~MultiResultStep() override {}

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return {}; }

  /// Add the buffer to the vector of kept buffers.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish does not do anything.
  void finish() override { getNextStep()->finish(); }

  /// Show the step parameters.
  /// It does nothing.
  void show(std::ostream&) const override{};

  /// Get the result.
  const std::vector<std::unique_ptr<base::DPBuffer>>& get() const {
    return buffers_;
  }
  std::vector<std::unique_ptr<base::DPBuffer>>& get() { return buffers_; }

  /// Get the size of the result.
  size_t size() const { return size_; }

  /// Clear the buffers.
  void clear() { size_ = 0; }

 private:
  std::vector<std::unique_ptr<base::DPBuffer>> buffers_;
  size_t size_;
};

}  // namespace steps
}  // namespace dp3

#endif
