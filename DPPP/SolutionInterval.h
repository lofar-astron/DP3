// SolutionInterval.h
//
// Copyright (C) 2020
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
/// @brief Solution Interval that can buffer multiple DPBuffers and data that is
/// relevant for using it in DDECal.
/// @author Lars Krombeen

#ifndef COMMON_SOLUTION_INTERVAL
#define COMMON_SOLUTION_INTERVAL

#include "../DPPP/DPBuffer.h"
#include "../DPPP/DPInput.h"

#include <casacore/casa/Arrays/Cube.h>

namespace DP3 {

namespace DPPP {

class SolutionInterval {
 public:
  SolutionInterval(DPInput* input, const std::size_t n_solution,
                   const std::size_t buffer_size, const std::size_t n_dir,
                   NSTimer timer);
  ~SolutionInterval();

  // Copy a buffer to the Solution Interval.
  void CopyBuffer(const DPBuffer&);

  /// Restore the flags and weights of added buffers to the original values
  void RestoreFlagsAndWeights();

  /// Resizes the pointer arrays so they only contain filled data.
  void Fit();

  /// Return the number of added buffers.
  const std::size_t Size() const { return buffer_index_; };

  /**
   * \defgroup Getters
   */
  /**@{*/
  const std::size_t NSolution() const { return n_solution_; }
  std::vector<casacore::Cube<casacore::Complex>>& ModelData() {
    return model_data_;
  }
  std::vector<std::vector<casacore::Complex*>>& ModelDataPtrs() {
    return model_data_ptrs_;
  };
  std::vector<float*>& WeightPtrs() { return weight_ptrs_; }
  std::vector<casacore::Complex*>& DataPtrs() { return data_ptrs_; }
  std::vector<std::vector<std::vector<casacore::Complex>>>& IDGBuffers() {
    return idg_buffers_;
  }
  /**@}*/

  /**
   * \defgroup Operators
   */
  /**@{*/
  /// Get a DPBuffer at index @param i
  const DPBuffer& operator[](int i) const { return buffers_[i]; }
  /**@}*/

 private:
  const std::size_t buffer_size_;
  const std::size_t n_solution_;

  NSTimer timer_;   ///< Timer from the step that is using it for metrics
  DPInput* input_;  ///< Input of DP3
  std::size_t buffer_index_;  ///< Current index where to insert the next buffer
  std::vector<DPBuffer> buffers_;  ///< Vector of DPBuffer copies

  std::vector<casacore::Complex*> data_ptrs_;
  std::vector<float*> weight_ptrs_;
  std::vector<casacore::Cube<bool>> original_flags_;
  std::vector<casacore::Cube<float>> original_weights_;

  /// For each timeslot, a vector of nDir buffers, each of size nbl x nch x npol
  std::vector<std::vector<casacore::Complex*>> model_data_ptrs_;
  std::vector<casacore::Cube<casacore::Complex>> model_data_;

  std::vector<std::vector<std::vector<casacore::Complex>>> idg_buffers_;
};
}  // namespace DPPP
}  // namespace DP3

#endif
