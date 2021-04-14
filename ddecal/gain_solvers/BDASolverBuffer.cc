// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BDASolverBuffer.h"

#include "../../base/BDABuffer.h"

#include <cassert>

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace base {

void BDASolverBuffer::AssignAndWeight(
    const std::vector<BDABuffer>& data_buffers,
    std::vector<std::vector<BDABuffer>>&& model_buffers) {
  const size_t n_buffers = model_buffers.front().size();
  data_.reserve(n_buffers);
  model_data_ = std::move(model_buffers);

  BDABuffer::Fields bda_fields;
  bda_fields.data_ = true;
  bda_fields.flags_ = false;
  bda_fields.weights_ = false;
  bda_fields.full_res_flags_ = false;

  for (size_t buffer = 0; buffer < n_buffers; ++buffer) {
    // Copy all unweighted data to data_.
    data_.emplace_back(data_buffers[buffer], bda_fields);
    BDABuffer& data_buffer = data_.back();

    const size_t n_rows = data_buffer.GetRows().size();
    for (size_t row = 0; row < n_rows; ++row) {
      const BDABuffer::Row& data_row = data_buffer.GetRows()[row];
      assert(kNCorrelations == data_row.n_correlations);

      for (size_t ch = 0; ch < data_row.n_channels; ++ch) {
        bool is_flagged = false;
        const size_t index = ch * kNCorrelations;
        const float* weights_ptr = data_buffers[buffer].GetWeights(row) + index;
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
        for (std::vector<BDABuffer>& direction_buffers : model_data_) {
          BDABuffer& model_buffer = direction_buffers[buffer];
          std::complex<float>* model_ptr = model_buffer.GetData(row) + index;
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
          for (std::vector<BDABuffer>& direction_buffers : model_data_) {
            BDABuffer& model_buffer = direction_buffers[buffer];
            std::complex<float>* model_ptr = model_buffer.GetData(row) + index;
            for (size_t cr = 0; cr < kNCorrelations; ++cr) {
              model_ptr[cr] = 0.0;
            }
          }
        }
      }
    }
  }
}

}  // namespace base
}  // namespace dp3
