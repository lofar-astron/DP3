// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaSolverBuffer.h"

#include <aocommon/matrix2x2.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#ifdef __SSE__
#include <immintrin.h>
#endif

using dp3::base::BdaBuffer;

namespace {
const size_t kNCorrelations = 4;

bool IsFinite(std::complex<float> c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}
}  // namespace

namespace dp3 {
namespace ddecal {

void BdaSolverBuffer::AppendAndWeight(
    std::unique_ptr<BdaBuffer> unweighted_buffer,
    const std::vector<std::string>& direction_names,
    bool keep_unweighted_model_data) {
  // Create a weighted_buffer with the main data and the model data for all
  // directions. Providing the 'data' field to the constructor would also work.
  // However, the constructor then copies all data buffers and not only those
  // in 'direction_names', resulting in more memory usage than needed.
  auto weighted_buffer =
      std::make_unique<BdaBuffer>(*unweighted_buffer, common::Fields());

  // Always copy the main data buffer.
  weighted_buffer->AddData();
  std::copy_n(unweighted_buffer->GetData(),
              unweighted_buffer->GetNumberOfElements(),
              weighted_buffer->GetData());

  // Copy or move the model buffer data into weighted_buffer.
  for (const std::string& name : direction_names) {
    if (keep_unweighted_model_data) {
      weighted_buffer->AddData(name);
      std::copy_n(unweighted_buffer->GetData(name),
                  unweighted_buffer->GetNumberOfElements(),
                  weighted_buffer->GetData(name));
    } else {
      weighted_buffer->MoveData(*unweighted_buffer, name, name);
    }
  }

  const size_t n_rows = weighted_buffer->GetRows().size();

  for (size_t row_index = 0; row_index < n_rows; ++row_index) {
    const BdaBuffer::Row& weighted_row = weighted_buffer->GetRows()[row_index];

    assert(weighted_row.interval <= time_interval_);
    assert(kNCorrelations == weighted_row.n_correlations);

    std::complex<float>* const row_unweighted_data =
        unweighted_buffer->GetData(row_index);
    std::complex<float>* const row_weighted_data =
        weighted_buffer->GetData(row_index);
    const float* const row_weights = unweighted_buffer->GetWeights(row_index);
    const bool* const row_flags = unweighted_buffer->GetFlags(row_index);

    std::vector<std::complex<float>*> unweighted_model_data;
    if (keep_unweighted_model_data) {
      unweighted_model_data.reserve(direction_names.size());
      for (const std::string& name : direction_names) {
        unweighted_model_data.push_back(
            unweighted_buffer->GetData(row_index, name));
      }
    }

    std::vector<std::complex<float>*> weighted_model_data;
    weighted_model_data.reserve(direction_names.size());
    for (const std::string& name : direction_names) {
      weighted_model_data.push_back(weighted_buffer->GetData(row_index, name));
    }

    for (size_t ch = 0; ch < weighted_row.n_channels; ++ch) {
      bool is_flagged = false;
      const size_t index = ch * kNCorrelations;
      std::complex<float>* const data_ptr = row_weighted_data + index;
      const bool* flag_ptr = row_flags + index;

      // Compute the square roots of the weights for the 4 correlations.
#ifdef __SSE__
      // All x86_64 CPUs support this SSE code, which computes the 4 square
      // roots in one instruction.
      const __m128 weights_sqrt =
          _mm_sqrt_ps(_mm_loadu_ps(row_weights + index));
#else
      // In case somebody want to build DP3 on ARM, this will work:
      const std::array<float, kNCorrelations> weights_sqrt{
          std::sqrt(row_weights[index + 0]), std::sqrt(row_weights[index + 1]),
          std::sqrt(row_weights[index + 2]), std::sqrt(row_weights[index + 3])};
#endif

      // Weight the 2x2 data matrix.
      for (size_t cr = 0; cr < kNCorrelations; ++cr) {
        is_flagged = is_flagged || !IsFinite(data_ptr[cr]) || flag_ptr[cr];
        data_ptr[cr] *= weights_sqrt[cr];
      }

      // Weight the model data.
      for (std::complex<float>* model_ptr : weighted_model_data) {
        model_ptr += index;

        for (size_t cr = 0; cr < kNCorrelations; ++cr) {
          is_flagged = is_flagged || !IsFinite(model_ptr[cr]);
          model_ptr[cr] *= weights_sqrt[cr];
        }
      }

      // If either the data or model data has non-finite values, set both the
      // data and model data to zero.
      if (is_flagged) {
        std::fill_n(data_ptr, kNCorrelations, 0.0);
        for (std::complex<float>* model_ptr : weighted_model_data) {
          std::fill_n(model_ptr + index, kNCorrelations, 0.0);
        }
      }
    }

    // Add row pointers to the correct solution interval.
    const int queue_index = RelativeIndex(weighted_row.time);
    assert(queue_index >= 0);

    // Add new solution intervals if needed.
    while (size_t(queue_index) >= data_rows_.Size()) {
      data_rows_.PushBack(std::vector<IntervalRow>());
    }

    // Add the row to the correct solution interval.
    std::vector<const std::complex<float>*> weighted_model_const_data(
        weighted_model_data.begin(), weighted_model_data.end());
    data_rows_[queue_index].emplace_back(
        weighted_row.time, weighted_row.baseline_nr, weighted_row.n_channels,
        weighted_row.n_correlations, row_unweighted_data, row_weighted_data,
        row_flags, row_weights, std::move(unweighted_model_data),
        std::move(weighted_model_const_data));

    int current_start_interval =
        RelativeIndex(weighted_row.time - weighted_row.interval / 2);
    assert(weighted_row.baseline_nr <
           last_complete_interval_per_baseline_.size());
    if (current_start_interval <
        last_complete_interval_per_baseline_[weighted_row.baseline_nr]) {
      throw std::invalid_argument(
          "Invalid input MS: BDA rows are not properly ordered.");
    }
    last_complete_interval_per_baseline_[weighted_row.baseline_nr] =
        current_start_interval - 1;
  }

  data_.PushBack(
      InputData{std::move(unweighted_buffer), std::move(weighted_buffer)});
}

void BdaSolverBuffer::Clear() {
  assert(!data_rows_.Empty());

  data_.Clear();
  data_rows_.Clear();
  done_.clear();

  // Add empty row vectors for the current solution interval.
  data_rows_.PushBack(std::vector<IntervalRow>());
}

void BdaSolverBuffer::AdvanceInterval() {
  assert(!data_rows_.Empty());

  data_rows_.PopFront();
  if (data_rows_.Empty()) data_rows_.PushBack(std::vector<IntervalRow>());

  ++current_interval_;
  for (int& bl_interval : last_complete_interval_per_baseline_) --bl_interval;

  // Remove old BdaBuffers.
  while (!data_.Empty()) {
    const bool all_rows_are_old =
        std::all_of(data_[0].unweighted->GetRows().begin(),
                    data_[0].unweighted->GetRows().end(),
                    [this](const BdaBuffer::Row& row) {
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

void BdaSolverBuffer::ApplySolutions(
    const std::vector<std::vector<std::complex<double>>>& solutions,
    const std::vector<double>& chan_block_start_freqs, size_t n_polarizations,
    const std::vector<int>& antennas1, const std::vector<int>& antennas2,
    const std::vector<std::vector<double>>& chan_freqs) {
  assert(!data_rows_.Empty());

  for (IntervalRow& row : data_rows_[0]) {
    // Map each (averaged) channel to a channel block.
    assert(!row.unweighted_model_data.empty());
    assert(chan_freqs[row.baseline_nr].size() == row.n_channels);
    std::vector<size_t> channel_blocks;
    channel_blocks.reserve(row.n_channels);
    size_t block = 0;
    for (const double freq : chan_freqs[row.baseline_nr]) {
      assert(freq > chan_block_start_freqs[block] &&
             freq < chan_block_start_freqs.back());
      while (chan_block_start_freqs[block + 1] < freq) ++block;
      channel_blocks.push_back(block);
    }

    const size_t n_directions = row.unweighted_model_data.size();
    const size_t ant1_index = antennas1[row.baseline_nr] * n_directions;
    const size_t ant2_index = antennas2[row.baseline_nr] * n_directions;

    for (size_t dir = 0; dir < n_directions; ++dir) {
      const size_t sol1_index = ant1_index + dir;
      const size_t sol2_index = ant2_index + dir;

      std::complex<float>* model_data_ptr = row.unweighted_model_data[dir];
      for (size_t ch = 0; ch < row.n_channels; ++ch) {
        const std::vector<std::complex<double>>& sol_block =
            solutions[channel_blocks[ch]];
        assert((sol1_index + 1) * n_polarizations <= sol_block.size());
        assert((sol2_index + 1) * n_polarizations <= sol_block.size());

        aocommon::MC2x2 model_data(model_data_ptr);
        switch (n_polarizations) {
          case 4: {
            const aocommon::MC2x2 sol1(&sol_block[sol1_index * 4]);
            const aocommon::MC2x2 sol2(&sol_block[sol2_index * 4]);
            model_data = (sol1 * model_data) * sol2.HermTranspose();
            break;
          }
          case 2: {
            const aocommon::MC2x2Diag sol1(&sol_block[sol1_index * 2]);
            const aocommon::MC2x2Diag sol2(&sol_block[sol2_index * 2]);
            model_data = (sol1 * model_data) * sol2.HermTranspose();
            break;
          }
          case 1: {
            const std::complex<double> sol_factor =
                sol_block[sol1_index] * std::conj(sol_block[sol2_index]);
            model_data = model_data * sol_factor;
            break;
          }
          default:
            assert(false);
        }
        model_data.AssignTo(model_data_ptr);
        model_data_ptr += kNCorrelations;
      }
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
