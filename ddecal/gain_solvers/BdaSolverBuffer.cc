// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaSolverBuffer.h"

#include <aocommon/matrix2x2.h>

#include <boost/make_unique.hpp>

#include <algorithm>
#include <cassert>
#include <limits>

using dp3::base::BDABuffer;

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace ddecal {

void BdaSolverBuffer::AppendAndWeight(
    std::unique_ptr<BDABuffer> unweighted_buffer,
    std::vector<std::unique_ptr<BDABuffer>>&& model_buffers) {
  const size_t n_directions = model_buffers.size();

  if (!data_.Empty() && n_directions != data_[0].model.size()) {
    throw std::invalid_argument("Model directions count does not match");
  }

  BDABuffer::Fields bda_fields(false);
  bda_fields.data = true;

  auto weighted_buffer =
      boost::make_unique<BDABuffer>(*unweighted_buffer, bda_fields);

  const size_t n_rows = weighted_buffer->GetRows().size();

  // Maximum start time of the row intervals.
  double max_start = std::numeric_limits<double>::min();

  for (size_t row = 0; row < n_rows; ++row) {
    const BDABuffer::Row& weighted_row = weighted_buffer->GetRows()[row];
    assert(kNCorrelations == weighted_row.n_correlations);

    for (size_t ch = 0; ch < weighted_row.n_channels; ++ch) {
      bool is_flagged = false;
      const size_t index = ch * kNCorrelations;
      const float* weights_ptr = unweighted_buffer->GetWeights(row) + index;
      const std::array<float, kNCorrelations> w_sqrt{
          std::sqrt(weights_ptr[0]), std::sqrt(weights_ptr[1]),
          std::sqrt(weights_ptr[2]), std::sqrt(weights_ptr[3])};

      // Weigh the 2x2 data matrix.
      std::complex<float>* data_ptr = weighted_row.data + index;
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
    const int queue_index = RelativeIndex(weighted_row.time);
    assert(queue_index >= 0);

    // Add new solution intervals if needed.
    while (size_t(queue_index) >= data_rows_.Size()) {
      AddInterval(data_rows_[0].model.size());
    }

    BDABuffer::Row& unweighted_row = unweighted_buffer->GetRows()[row];
    data_rows_[queue_index].unweighted.push_back(&unweighted_row);
    data_rows_[queue_index].weighted.push_back(&weighted_row);
    for (size_t dir = 0; dir < n_directions; ++dir) {
      data_rows_[queue_index].model[dir].push_back(
          &model_buffers[dir]->GetRows()[row]);
    }

    max_start =
        std::max(max_start, weighted_row.time - weighted_row.interval / 2);
  }

  data_.PushBack(InputData{std::move(unweighted_buffer),
                           std::move(weighted_buffer),
                           std::move(model_buffers)});

  // Update last_complete_interval_.
  int max_start_interval = RelativeIndex(max_start);
  assert(max_start_interval > last_complete_interval_);
  last_complete_interval_ = max_start_interval - 1;
}

void BdaSolverBuffer::Clear() {
  assert(!data_rows_.Empty());
  const size_t n_directions = data_rows_[0].model.size();

  data_.Clear();
  data_rows_.Clear();
  done_.clear();

  // Add empty row vectors for the current solution interval.
  AddInterval(n_directions);
}

void BdaSolverBuffer::AdvanceInterval() {
  assert(!data_rows_.Empty());
  const size_t n_directions = data_rows_[0].model.size();

  data_rows_.PopFront();
  if (data_rows_.Empty()) AddInterval(n_directions);

  ++current_interval_;
  --last_complete_interval_;

  // Remove old BDABuffers.
  while (!data_.Empty()) {
    const bool all_rows_are_old =
        std::all_of(data_[0].unweighted->GetRows().begin(),
                    data_[0].unweighted->GetRows().end(),
                    [this](const BDABuffer::Row& row) {
                      return RelativeIndex(row.time) < 0;
                    });

    if (all_rows_are_old) {
      done_.push_back(std::move(data_[0].unweighted));
      data_.PopFront();
    } else {
      break;
    }
  }
}

void BdaSolverBuffer::AddInterval(size_t n_directions) {
  data_rows_.PushBack(
      {std::vector<BDABuffer::Row*>(), std::vector<const BDABuffer::Row*>(),
       std::vector<std::vector<const base::BDABuffer::Row*>>(n_directions)});
}

void BdaSolverBuffer::SubtractCorrectedModel(
    const std::vector<std::vector<std::complex<float>>>& solutions,
    size_t n_channel_blocks, bool full_jones, const std::vector<int>& antennas1,
    const std::vector<int>& antennas2) {
  // data_ and data_rows_ still hold the original unweighted input data, since
  // the Solver doesn't change those. Here we apply the solutions to all the
  // model data directions and subtract them from the unweighted input data.
  assert(!data_rows_.Empty());
  const size_t n_directions = data_rows_[0].model.size();

  for (size_t row = 0; row < data_rows_[0].unweighted.size(); ++row) {
    BDABuffer::Row* unweighted_row = data_rows_[0].unweighted[row];

    // Determine channel block endings for this row.
    std::vector<size_t> chan_block_end;
    chan_block_end.reserve(n_channel_blocks);
    for (size_t b = 0; b < n_channel_blocks; ++b) {
      chan_block_end.push_back((b + 1) * unweighted_row->n_channels /
                               n_channel_blocks);
    }

    const size_t ant1_index =
        antennas1[unweighted_row->baseline_nr] * n_directions;
    const size_t ant2_index =
        antennas2[unweighted_row->baseline_nr] * n_directions;

    for (size_t dir = 0; dir < n_directions; ++dir) {
      const BDABuffer::Row* model_row = data_rows_[0].model[dir][row];
      assert(model_row->n_correlations == kNCorrelations);

      const size_t sol1_index = ant1_index + dir;
      const size_t sol2_index = ant2_index + dir;

      size_t chan_block = 0;
      std::complex<float>* unweighted_data = unweighted_row->data;
      const std::complex<float>* model_data = model_row->data;
      for (size_t ch = 0; ch < unweighted_row->n_channels; ++ch) {
        if (ch == chan_block_end[chan_block]) ++chan_block;
        const std::vector<std::complex<float>>& sol_block =
            solutions[chan_block];

        if (full_jones) {
          assert((sol1_index + 1) * kNCorrelations <= sol_block.size());
          assert((sol2_index + 1) * kNCorrelations <= sol_block.size());
          const aocommon::MC2x2 sol1(&sol_block[sol1_index * kNCorrelations]);
          const aocommon::MC2x2 sol2(&sol_block[sol2_index * kNCorrelations]);
          const aocommon::MC2x2 solved =
              sol1.Multiply(aocommon::MC2x2(model_data)).MultiplyHerm(sol2);
          for (size_t corr = 0; corr < kNCorrelations; ++corr) {
            unweighted_data[corr] -= solved[corr];
          }
        } else {
          assert(sol1_index < sol_block.size());
          assert(sol2_index < sol_block.size());
          const std::complex<float> sol_factor =
              sol_block[sol1_index] * std::conj(sol_block[sol2_index]);
          for (size_t corr = 0; corr < kNCorrelations; ++corr) {
            unweighted_data[corr] -= sol_factor * model_data[corr];
          }
        }

        unweighted_data += kNCorrelations;
        model_data += kNCorrelations;
      }
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
