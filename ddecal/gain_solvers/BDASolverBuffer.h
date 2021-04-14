// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_SOLVER_BUFFER_H
#define DDECAL_BDA_SOLVER_BUFFER_H

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class BDABuffer;

class BDASolverBuffer {
 public:
  BDASolverBuffer() : data_(), model_data_() {}

  /**
   * This function takes (unweighted) data and model data, as well as a
   * weights array, weights these and initializes the buffer from
   * them.
   *
   * @param data_buffers A vector with buffers that cover the solution interval.
   * Each buffer should contain unweighted data and the corresponding weight
   * values.
   * @param model_buffers Buffers with model data. For each direction,
   * model_buffers holds a vector with buffers with the unweighted model
   * data. The BDA layout of the buffer in the inner vector should match the
   * layout of 'data_buffers'. The BDASolverBuffer takes ownership of the model
   * buffers.
   */
  void AssignAndWeight(const std::vector<BDABuffer>& data_buffers,
                       std::vector<std::vector<BDABuffer>>&& model_buffers);

  /**
   * @return The weighted data buffers.
   */
  const std::vector<BDABuffer>& GetData() const { return data_; }

  /**
   * @param direction Direction index.
   * @return The weighted model buffers for the given direction.
   */
  const std::vector<BDABuffer>& GetModelData(size_t direction) const {
    return model_data_[direction];
  }

 private:
  std::vector<BDABuffer> data_;
  std::vector<std::vector<BDABuffer>> model_data_;
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_BDA_SOLVER_BUFFER_H
