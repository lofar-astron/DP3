// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_TEST_UNIT_MOCKSTEP_H_
#define DP3_STEPS_TEST_UNIT_MOCKSTEP_H_

#include <functional>
#include <vector>

#include "ThrowStep.h"

namespace dp3 {
namespace steps {

class MockStep : public test::ThrowStep {
 public:
  MockStep() : bda_buffers_(), regular_buffers_(), finish_count_(0) {}

  ~MockStep() override {}

  /**
   * Mocked process() function for regular buffers.
   * Adds the regular buffer to an internal list. Use GetRegularBuffers for
   * accessing them.
   */
  bool process(std::unique_ptr<base::DPBuffer> buffer) override {
    regular_buffers_.push_back(std::move(buffer));
    return true;
  }

  /**
   * Mocked process() function for bda buffers.
   * Adds the bda buffer to an internal list. Use GetBdaBuffers for
   * accessing them.
   */
  bool process(std::unique_ptr<base::BdaBuffer> buffer) override {
    bda_buffers_.push_back(std::move(buffer));
    return true;
  }

  /**
   * Mocked finish() function, which counts the number of calls.
   * Use FinishCount() for accessing the count.
   */
  void finish() override { ++finish_count_; }

  const std::vector<std::unique_ptr<base::BdaBuffer>>& GetBdaBuffers() const {
    return bda_buffers_;
  }

  const std::vector<std::unique_ptr<base::DPBuffer>>& GetRegularBuffers()
      const {
    return regular_buffers_;
  }

  void ClearBdaBuffers() { bda_buffers_.clear(); }

  std::size_t FinishCount() const { return finish_count_; };

  std::size_t TotalRowCount() const;

 private:
  std::vector<std::unique_ptr<base::BdaBuffer>> bda_buffers_;
  std::vector<std::unique_ptr<base::DPBuffer>> regular_buffers_;
  std::size_t finish_count_;
};

}  // namespace steps
}  // namespace dp3

#endif
