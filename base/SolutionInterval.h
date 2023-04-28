// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @brief Solution Interval that can buffer multiple DPBuffers and data that is
/// relevant for using it in DDECal.
/// @author Lars Krombeen

#ifndef DP3_BASE_SOLUTION_INTERVAL_
#define DP3_BASE_SOLUTION_INTERVAL_

#include <casacore/casa/Arrays/Cube.h>

#include <dp3/base/DPBuffer.h>

namespace dp3 {
namespace base {

class SolutionInterval {
 public:
  explicit SolutionInterval(std::size_t buffer_size);

  // Append a buffer to the Solution Interval.
  void PushBack(std::unique_ptr<DPBuffer> buffer);

  /// Return the number of added buffers.
  std::size_t Size() const { return buffers_.size(); };

  /**
   * \defgroup Getters
   */
  /**@{*/
  std::vector<std::unique_ptr<DPBuffer>>& DataBuffers() { return buffers_; }
  /**@}*/

  /**
   * \defgroup Operators
   */
  /**@{*/
  /// Get a DPBuffer at index @param i
  const DPBuffer& operator[](int i) const { return *buffers_[i]; }
  /**@}*/

 private:
  std::vector<std::unique_ptr<DPBuffer>> buffers_;
};

}  // namespace base
}  // namespace dp3

#endif
