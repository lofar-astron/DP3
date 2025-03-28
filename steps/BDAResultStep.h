// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BDARESULTSTEP_H
#define DP3_BDARESULTSTEP_H

#include <dp3/steps/Step.h>
#include <dp3/base/BdaBuffer.h>

#include <cassert>
#include <vector>

namespace dp3 {
namespace steps {

/**
 * Defines a result step for BDA buffers.
 * - The process() function adds all buffers to an internal queue.
 * - The result step should always be the last step in a series of steps.
 * - The owner of the BDAResultStep (BDADDECal) extracts buffers from the queue.
 */
class BDAResultStep : public Step {
 public:
  /// Creates an empty BDAResultStep.
  BDAResultStep() : buffers_() {}

  ~BDAResultStep() override {}

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return {}; }

  /// Adds a buffer to the internal queue.
  bool process(std::unique_ptr<base::BdaBuffer> buffer) override {
    assert(!getNextStep());  // The result step should be the last step.
    buffers_.push_back(std::move(buffer));
    return true;
  }

  /// Does nothing.
  void finish() override {}

  /// Does nothing.
  void show(std::ostream&) const override {}

  /// Extracts all stored buffers from the result step.
  std::vector<std::unique_ptr<base::BdaBuffer>> Extract() {
    std::vector<std::unique_ptr<base::BdaBuffer>> result;
    result.swap(buffers_);
    return result;
  }

 private:
  std::vector<std::unique_ptr<base::BdaBuffer>> buffers_;
};

}  // namespace steps
}  // namespace dp3

#endif
