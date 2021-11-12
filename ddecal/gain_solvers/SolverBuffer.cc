// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBuffer.h"

#include "../../base/DPBuffer.h"

#include <cassert>

using dp3::base::DPBuffer;

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace ddecal {

SolverBuffer::SolverBuffer()
    : n_baselines_(0), n_channels_(0), data_(), model_buffers_() {}

void SolverBuffer::AssignAndWeight(
    const std::vector<DPBuffer>& unweighted_data_buffers,
    std::vector<std::vector<DPBuffer>>&& model_buffers) {
  const size_t n_times = model_buffers.size();
  assert(unweighted_data_buffers.size() >= n_times);
  data_.resize(n_times);
  model_buffers_ = std::move(model_buffers);

  for (size_t timestep = 0; timestep != n_times; ++timestep) {
    const casacore::Cube<std::complex<float>>& unweighted_data =
        unweighted_data_buffers[timestep].getData();
    const casacore::Cube<float>& weights =
        unweighted_data_buffers[timestep].getWeights();
    assert(unweighted_data.shape() == weights.shape());
    assert(timestep == 0 ||
           unweighted_data.shape() ==
               unweighted_data_buffers.front().getData().shape());

    n_baselines_ = unweighted_data.shape()[2];
    n_channels_ = unweighted_data.shape()[1];
    assert(kNCorrelations == unweighted_data.shape()[0]);

    data_[timestep].resize(unweighted_data.size());
    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      for (size_t ch = 0; ch < n_channels_; ++ch) {
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
        for (DPBuffer& model_buffer : model_buffers_[timestep]) {
          Complex* model_ptr = model_buffer.getData().data();
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
          for (DPBuffer& model_buffer : model_buffers_[timestep]) {
            Complex* model_ptr = model_buffer.getData().data();
            for (size_t cr = 0; cr < kNCorrelations; ++cr) {
              model_ptr[index + cr] = 0.0;
            }
          }
        }
      }
    }
  }
}

const std::complex<float>* SolverBuffer::DataPointer(size_t time_index,
                                                     size_t baseline,
                                                     size_t channel) const {
  return &data_[time_index]
               [(baseline * n_channels_ + channel) * kNCorrelations];
}

const std::complex<float>* SolverBuffer::ModelDataPointer(
    size_t time_index, size_t direction, size_t baseline,
    size_t channel) const {
  const DPBuffer& buffer = model_buffers_[time_index][direction];
  return &buffer.getData()(0, channel, baseline);
}

void SolverBuffer::CopyDataChannels(size_t time_index, size_t channel_begin,
                                    size_t channel_end,
                                    std::complex<float>* destination) const {
  const size_t baseline_size = n_channels_ * kNCorrelations;
  const size_t channels_size = (channel_end - channel_begin) * kNCorrelations;

  const std::complex<float>* data = DataPointer(time_index, 0, channel_begin);
  for (size_t bl = 0; bl < n_baselines_; ++bl) {
    std::copy_n(data, channels_size, destination);
    data += baseline_size;
    destination += channels_size;
  }
}

}  // namespace ddecal
}  // namespace dp3
