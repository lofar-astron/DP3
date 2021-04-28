// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BDASolverBuffer.h"

#include <cassert>

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace base {

void BDASolverBuffer::AppendAndWeight(
    const BDABuffer& input_data_buffer,
    std::vector<std::unique_ptr<BDABuffer>>&& model_buffers) {
  if (time_interval_ == 0.0)
    throw std::logic_error("Solution interval is not set.");

  if (model_data_.size() != model_buffers.size())
    throw std::invalid_argument("Invalid number of model directions.");

  BDABuffer::Fields bda_fields(false);
  bda_fields.data = true;

  // Copy all unweighted data to data_.
  data_.emplace_back(input_data_buffer, bda_fields);
  BDABuffer& data_buffer = data_.back();

  const size_t n_rows = data_buffer.GetRows().size();
  for (size_t row = 0; row < n_rows; ++row) {
    const BDABuffer::Row& data_row = data_buffer.GetRows()[row];
    assert(kNCorrelations == data_row.n_correlations);

    for (size_t ch = 0; ch < data_row.n_channels; ++ch) {
      bool is_flagged = false;
      const size_t index = ch * kNCorrelations;
      const float* weights_ptr = input_data_buffer.GetWeights(row) + index;
      const std::array<float, kNCorrelations> w_sqrt{
          std::sqrt(weights_ptr[0]), std::sqrt(weights_ptr[1]),
          std::sqrt(weights_ptr[2]), std::sqrt(weights_ptr[3])};

      // Weigh the 2x2 data matrix.
      std::complex<float>* data_ptr = data_row.data + index;
      for (size_t cr = 0; cr < kNCorrelations; ++cr) {
        is_flagged = is_flagged || !IsFinite(data_ptr[cr]);
        data_ptr[cr] *= w_sqrt[cr];
      }

      // Weigh the model data.
      for (std::unique_ptr<BDABuffer>& model_buffer : model_buffers) {
        std::complex<float>* model_ptr = model_buffer->GetData(row) + index;
        for (size_t cr = 0; cr < kNCorrelations; ++cr) {
          is_flagged = is_flagged || !IsFinite(model_ptr[cr]);
          model_ptr[cr] *= w_sqrt[cr];
        }
      }

      // If either the data or model data has non-finite values, set both the
      // data and model data to zero.
      if (is_flagged) {
        for (size_t cr = 0; cr < kNCorrelations; ++cr) {
          data_ptr[cr] = 0.0;
        }
        for (std::unique_ptr<BDABuffer>& model_buffer : model_buffers) {
          std::complex<float>* model_ptr = model_buffer->GetData(row) + index;
          for (size_t cr = 0; cr < kNCorrelations; ++cr) {
            model_ptr[cr] = 0.0;
          }
        }
      }
    }
  }

  // Append the model_buffers to the list for each direction.
  const size_t n_directions = model_data_.size();
  for (size_t dir = 0; dir < n_directions; ++dir) {
    model_data_[dir].push_back(std::move(model_buffers[dir]));
  }
}

}  // namespace base
}  // namespace dp3
