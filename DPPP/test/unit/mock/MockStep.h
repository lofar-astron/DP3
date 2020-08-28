// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief DPStep mock, for use in tests.
/// @author Maik Nijhuis

#include "../../../DPStep.h"

#include <functional>
#include <vector>

namespace DP3 {
namespace DPPP {

class MockStep : public DPStep {
 public:
  /**
   * Constructor.
   * @param check_buffer The function that must be called when process(DPBuffer)
   *        is called. Use nullptr when the mock does not expect these calls.
   */
  explicit MockStep(
      std::function<void(const DPBuffer&)>* check_buffer = nullptr);

  ~MockStep() override;

  /**
   * Mocked process() function for regular buffers.
   * If no check_buffer function was set, the unit test fails.
   * Otherwise, this function forwards the buffer to the check_buffer function.
   */
  bool process(const DPBuffer&) override;

  /**
   * Mocked process() function for bda buffers.
   * Adds the bda buffer to an internal list. Use getBdaBuffers for
   * accessing them.
   */
  bool process(std::unique_ptr<BDABuffer>) override;

  /**
   * Mocked finish() function, which counts the number of calls.
   * Use getFinishCount() for accessing the count.
   */
  void finish() override { ++finish_count_; }

  /**
   * Mocked show() function. The unit test fails if it is called.
   */
  void show(std::ostream&) const override;

  const std::vector<std::unique_ptr<BDABuffer>>& GetBdaBuffers() const {
    return bda_buffers_;
  }

  void ClearBdaBuffers();

  std::size_t FinishCount() const { return finish_count_; };

  std::size_t TotalRowCount() const;

 private:
  std::function<void(const DPBuffer&)>* check_buffer_;
  std::vector<std::unique_ptr<BDABuffer>> bda_buffers_;
  std::size_t finish_count_;
};

}  // namespace DPPP
}  // namespace DP3
