// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_BUFFER_H
#define DDECAL_SOLVER_BUFFER_H

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class DPBuffer;

class SolverBuffer {
 public:
  typedef std::complex<float> Complex;

  explicit SolverBuffer(size_t n_directions = 0)
      : n_directions_(n_directions) {}

  void SetDirectionCount(size_t n_directions) { n_directions_ = n_directions; }

  /**
   * This function takes (unweighted) data and model data, as well as a
   * weights array, and weights these and initializes the buffer from
   * this.
   *
   * The data ordering is defined in SolverBase.
   */
  void AssignAndWeight(
      const std::vector<DPBuffer>& unweighted_data_buffers,
      const std::vector<std::vector<DPBuffer*>>& model_buffers);

  const std::vector<std::vector<Complex>>& Data() const { return data_; }

 private:
  size_t n_directions_;

  std::vector<std::vector<Complex>> data_;
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_SOLVER_BUFFER_H
