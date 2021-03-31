// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBuffer.h"

namespace dp3 {
namespace base {

void SolverBuffer::AssignAndWeight(
    const std::vector<Complex*>& unweighted_data,
    const std::vector<float*>& weights,
    std::vector<std::vector<Complex*>>&& unweighted_model_data) {
  const size_t n_times = unweighted_data.size();
  data_.resize(n_times);
  model_data_ = std::move(unweighted_model_data);

  for (size_t timestep = 0; timestep != n_times; ++timestep) {
    data_[timestep].resize(n_baselines_ * n_channels_ * 4);
    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      for (size_t ch = 0; ch != n_channels_; ++ch) {
        bool isFlagged = false;

        // Copy and weigh the 2x2 data matrix
        for (size_t cr = 0; cr < 4; ++cr) {
          const size_t index = (bl * n_channels_ + ch) * 4 + cr;

          if (!Isfinite(unweighted_data[timestep][index]))
            isFlagged = true;
          else {
            float wSqrt = std::sqrt(weights[timestep][index]);
            data_[timestep][index] = unweighted_data[timestep][index] * wSqrt;
          }
        }

        // Weigh the model data. This is done in a separate loop to loop
        // over the data contiguously in memory.
        const size_t index = (bl * n_channels_ + ch) * 4;
        for (size_t dir = 0; dir < n_directions_; ++dir) {
          for (size_t cr = 0; cr < 4; ++cr) {
            float wSqrt = std::sqrt(weights[timestep][index + cr]);

            if (!Isfinite(model_data_[timestep][dir][index + cr]))
              isFlagged = true;
            model_data_[timestep][dir][index + cr] *= wSqrt;
          }
        }

        // If either the data or model data has non-finite values, both the
        // data and model data are set to zero.
        if (isFlagged) {
          const size_t index = (bl * n_channels_ + ch) * 4;
          for (size_t cr = 0; cr < 4; ++cr) {
            data_[timestep][index + cr] = 0.0;
          }
          for (size_t dir = 0; dir < n_directions_; ++dir) {
            for (size_t cr = 0; cr < 4; ++cr) {
              model_data_[timestep][dir][index + cr] = 0.0;
            }
          }
        }
      }
    }
  }
}

}  // namespace base
}  // namespace dp3
