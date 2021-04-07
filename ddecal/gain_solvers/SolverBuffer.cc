// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBuffer.h"

#include "../../base/DPBuffer.h"

#include <cassert>

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace base {

void SolverBuffer::AssignAndWeight(
    const std::vector<DPBuffer>& unweighted_data_buffers,
    const std::vector<std::vector<DPBuffer*>>& model_buffers) {
  const size_t n_times = model_buffers.size();
  data_.resize(n_times);

  for (size_t timestep = 0; timestep != n_times; ++timestep) {
    const casacore::Cube<std::complex<float>>& unweighted_data =
        unweighted_data_buffers[timestep].getData();
    const casacore::Cube<float>& weights =
        unweighted_data_buffers[timestep].getWeights();
    const size_t n_baselines = unweighted_data.shape()[2];
    const size_t n_channels = unweighted_data.shape()[1];
    assert(kNCorrelations == unweighted_data.shape()[0]);

    data_[timestep].resize(unweighted_data.size());
    for (size_t bl = 0; bl < n_baselines; ++bl) {
      for (size_t ch = 0; ch < n_channels; ++ch) {
        const size_t index =
            &unweighted_data(0, ch, bl) - unweighted_data.data();
        bool is_flagged = false;
        const std::array<float, kNCorrelations> w_sqrt{
            std::sqrt(weights.data()[index + 0]),
            std::sqrt(weights.data()[index + 1]),
            std::sqrt(weights.data()[index + 2]),
            std::sqrt(weights.data()[index + 3])};

        // Copy and weigh the 2x2 data matrix
        for (size_t cr = 0; cr < kNCorrelations; ++cr) {
          const std::complex<float> ud = unweighted_data.data()[index + cr];
          is_flagged = is_flagged || !IsFinite(ud);
          data_[timestep][index + cr] = ud * w_sqrt[cr];
        }

        // Weigh the model data. This is done in a separate loop to loop
        // over the data contiguously in memory.
        for (size_t dir = 0; dir < n_directions_; ++dir) {
          Complex* model_ptr = model_buffers[timestep][dir]->getData().data();
          for (size_t cr = 0; cr < kNCorrelations; ++cr) {
            is_flagged = is_flagged || !IsFinite(model_ptr[index + cr]);
            model_ptr[index + cr] *= w_sqrt[cr];
          }
        }

        // If either the data or model data has non-finite values, set both the
        // data and model data to zero.
        if (is_flagged) {
          for (size_t cr = 0; cr < kNCorrelations; ++cr) {
            data_[timestep][index + cr] = 0.0;
          }
          for (size_t dir = 0; dir < n_directions_; ++dir) {
            Complex* model_ptr = model_buffers[timestep][dir]->getData().data();
            for (size_t cr = 0; cr < kNCorrelations; ++cr) {
              model_ptr[index + cr] = 0.0;
            }
          }
        }
      }
    }
  }
}

}  // namespace base
}  // namespace dp3
