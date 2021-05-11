// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step mock, for use in tests.
/// @author Maik Nijhuis

#ifndef MOCK_STEP_H
#define MOCK_STEP_H

#include "../../../Step.h"

#include <functional>
#include <vector>

namespace dp3 {
namespace steps {

class MockStep : public Step {
 public:
  /**
   * Constructor.
   * @param check_buffer The function that must be called when process(DPBuffer)
   *        is called. Use nullptr when the mock does not expect these calls.
   */
  explicit MockStep(
      std::function<void(const base::DPBuffer &)> *check_buffer = nullptr);

  ~MockStep() override;

  /**
   * Mocked process() function for regular buffers.
   * If no check_buffer function was set, the unit test fails.
   * Otherwise, this function forwards the buffer to the check_buffer function.
   */
  bool process(const base::DPBuffer &) override;

  /**
   * Mocked process() function for bda buffers.
   * Adds the bda buffer to an internal list. Use getBdaBuffers for
   * accessing them.
   */
  bool process(std::unique_ptr<base::BDABuffer>) override;

  /**
   * Mocked finish() function, which counts the number of calls.
   * Use getFinishCount() for accessing the count.
   */
  void finish() override { ++finish_count_; }

  /**
   * Mocked show() function. The unit test fails if it is called.
   */
  void show(std::ostream &) const override;

  const std::vector<std::unique_ptr<base::BDABuffer>> &GetBdaBuffers() const {
    return bda_buffers_;
  }

  const std::vector<base::DPBuffer> &GetRegularBuffers() const {
    return regular_buffers_;
  }

  void ClearBdaBuffers();

  std::size_t FinishCount() const { return finish_count_; };

  std::size_t TotalRowCount() const;

 private:
  std::function<void(const base::DPBuffer &)> *check_buffer_;
  std::vector<std::unique_ptr<base::BDABuffer>> bda_buffers_;
  std::vector<base::DPBuffer> regular_buffers_;
  std::size_t finish_count_;
};

}  // namespace steps
}  // namespace dp3

#endif
