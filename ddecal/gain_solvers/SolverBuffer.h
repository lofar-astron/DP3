// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_BUFFER_H
#define DDECAL_SOLVER_BUFFER_H

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class SolverBuffer {
 public:
  typedef std::complex<float> Complex;

  SolverBuffer() : n_directions_(0), n_channels_(0), n_baselines_(0) {}

  SolverBuffer(size_t n_directions, size_t n_channels, size_t n_baselines)
      : n_directions_(n_directions),
        n_channels_(n_channels),
        n_baselines_(n_baselines) {}

  void SetDimensions(size_t n_directions, size_t n_channels,
                     size_t n_baselines) {
    n_directions_ = n_directions;
    n_channels_ = n_channels;
    n_baselines_ = n_baselines;
  }

  /**
   * This function takes (unweighted) data and model data, as well as a
   * weights array, and weights these and initializes the buffer from
   * this.
   *
   * Note the confusing semantics: the unweighted_data and weights arrays
   * are copied (without changing the input) whereas unweighted_model_data
   * is moved from, and therefore left empty after the call. These
   * semantics help in preventing copies inside DDECal.
   *
   * The data ordering is defined in SolverBase.
   */
  void AssignAndWeight(
      const std::vector<Complex*>& unweighted_data,
      const std::vector<float*>& weights,
      std::vector<std::vector<Complex*>>&& unweighted_model_data);

  const std::vector<std::vector<Complex>>& Data() const { return data_; }

  const std::vector<std::vector<Complex*>>& ModelData() const {
    return model_data_;
  }

 private:
  size_t n_directions_, n_channels_, n_baselines_;

  std::vector<std::vector<Complex>> data_;
  std::vector<std::vector<Complex*>> model_data_;

  static bool Isfinite(Complex c) {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_SOLVER_BUFFER_H
