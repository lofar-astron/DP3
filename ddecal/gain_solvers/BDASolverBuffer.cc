// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BDASolverBuffer.h"

#include <boost/make_unique.hpp>

#include <cassert>
#include <limits>

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
  const size_t n_directions = model_data_.size();

  if (n_directions != model_buffers.size())
    throw std::invalid_argument("Invalid number of model directions.");

  BDABuffer::Fields bda_fields(false);
  bda_fields.data = true;

  // Copy all unweighted data to data_.
  BDABuffer& data_buffer = *data_.PushBack(
      boost::make_unique<BDABuffer>(input_data_buffer, bda_fields));
  const size_t n_rows = data_buffer.GetRows().size();

  // Maximum start time of the row intervals.
  double max_start = std::numeric_limits<double>::min();

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

    // Add row pointers to the correct solution interval.
    assert(data_row.time > time_start_);
    const size_t interval_index =
        (data_row.time - time_start_) / time_interval_;

    // Add new solution intervals if needed.
    while (interval_index >= data_rows_.Size()) AddInterval();

    data_rows_[interval_index].push_back(&data_row);
    for (size_t dir = 0; dir < n_directions; ++dir) {
      model_rows_[dir][interval_index].push_back(
          &model_buffers[dir]->GetRows()[row]);
    }

    max_start = std::max(max_start, data_row.time - data_row.interval / 2);
  }

  // Append the model_buffers to the list for each direction.
  for (size_t dir = 0; dir < n_directions; ++dir) {
    model_data_[dir].PushBack(std::move(model_buffers[dir]));
  }

  // Update last_complete_interval_.
  int max_start_interval = (max_start - time_start_) / time_interval_;
  assert(max_start_interval > last_complete_interval_);
  last_complete_interval_ = max_start_interval - 1;
}

void BDASolverBuffer::Clear() {
  data_.Clear();
  data_rows_.Clear();

  for (auto& model_buffer_queue : model_data_) model_buffer_queue.Clear();
  for (auto& model_row_queue : model_rows_) model_row_queue.Clear();

  // Add empty row vectors for the current solution interval.
  AddInterval();
}

void BDASolverBuffer::AdvanceInterval() {
  assert(!data_rows_.Empty());

  data_rows_.PopFront();
  for (auto& model_row_queue : model_rows_) model_row_queue.PopFront();

  if (data_rows_.Empty()) AddInterval();

  time_start_ += time_interval_;
  --last_complete_interval_;

  // Remove old BDABuffers.
  while (!data_.Empty()) {
    bool all_rows_are_old = true;
    for (const BDABuffer::Row& row : data_[0]->GetRows()) {
      if (BDABuffer::TimeIsGreaterEqual(row.time, time_start_)) {
        all_rows_are_old = false;
        break;
      }
    }
    if (all_rows_are_old) {
      data_.PopFront();
      for (auto& model_data_queue : model_data_) model_data_queue.PopFront();
    } else {
      break;
    }
  }
}

void BDASolverBuffer::AddInterval() {
  data_rows_.PushBack(std::vector<const BDABuffer::Row*>());
  for (auto& model_row_queue : model_rows_) {
    model_row_queue.PushBack(std::vector<const BDABuffer::Row*>());
  }
}

}  // namespace base
}  // namespace dp3
