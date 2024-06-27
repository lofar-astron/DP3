// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolveData.h"

#include "BdaSolverBuffer.h"

#include <xtensor/xview.hpp>

#include <cassert>

using dp3::base::BDABuffer;

namespace dp3 {
namespace ddecal {

template <typename MatrixType>
MatrixType ToSpecificMatrix(const std::complex<float>* data);
template <>
aocommon::MC2x2F ToSpecificMatrix(const std::complex<float>* data) {
  return aocommon::MC2x2F(data);
}
template <>
aocommon::MC2x2FDiag ToSpecificMatrix(const std::complex<float>* data) {
  return aocommon::MC2x2FDiag(data[0], data[3]);
}

template <typename MatrixType>
SolveData<MatrixType>::SolveData(
    const std::vector<base::DPBuffer>& buffers,
    const std::vector<std::string>& direction_names, size_t n_channel_blocks,
    size_t n_antennas, const std::vector<size_t>& n_solutions_per_direction,
    const std::vector<int>& antennas1, const std::vector<int>& antennas2)
    : channel_blocks_(n_channel_blocks) {
  std::vector<size_t> channel_begin(n_channel_blocks + 1, 0);

  const size_t n_times = buffers.size();
  const size_t n_baselines_in_buffers =
      buffers.empty() ? 0 : buffers.front().GetData().shape(0);
  const size_t n_channels =
      buffers.empty() ? 0 : buffers.front().GetData().shape(1);
  const size_t n_directions = direction_names.size();
  const bool has_weights = buffers.front().GetWeights().size() != 0;

  // Count nr of baselines with different antennas.
  size_t n_baselines = 0;
  for (size_t baseline = 0; baseline < n_baselines_in_buffers; ++baseline) {
    assert(size_t(antennas1[baseline]) < n_antennas &&
           size_t(antennas2[baseline]) < n_antennas);
    if (antennas1[baseline] != antennas2[baseline]) ++n_baselines;
  }

  // Initialize n_solutions_ of the first channel block as template
  // for the other channel blocks
  channel_blocks_.front().n_solutions_.reserve(n_directions);
  for (size_t direction = 0; direction != n_directions; ++direction) {
    channel_blocks_.front().n_solutions_.emplace_back(
        std::min(n_solutions_per_direction[direction], n_times));
  }

  // Count nr of visibilities per channel block and allocate memory.
  for (size_t channel_block_index = 0; channel_block_index != n_channel_blocks;
       ++channel_block_index) {
    ChannelBlockData& cb_data = channel_blocks_[channel_block_index];

    channel_begin[channel_block_index + 1] =
        (channel_block_index + 1) * n_channels / n_channel_blocks;
    const size_t channel_block_size = channel_begin[channel_block_index + 1] -
                                      channel_begin[channel_block_index];

    const size_t n_visibilities = n_times * n_baselines * channel_block_size;
    cb_data.Resize(n_visibilities, n_directions);
    if (has_weights) cb_data.ResizeWeights(n_visibilities);

    cb_data.n_solutions_ = channel_blocks_.front().n_solutions_;
  }

  // Construct an array that identifies the start solution index per direction
  std::vector<size_t> solution_start_indices;
  solution_start_indices.reserve(n_directions);
  size_t solution_start_counter = 0;
  for (size_t direction = 0; direction != n_directions; ++direction) {
    solution_start_indices.emplace_back(solution_start_counter);
    solution_start_counter += channel_blocks_.front().n_solutions_[direction];
  }

  // Fill all channel blocks with data.
  std::vector<size_t> visibility_indices(n_channel_blocks, 0);
  for (size_t time_index = 0; time_index < n_times; ++time_index) {
    const base::DPBuffer::DataType& data = buffers[time_index].GetData("");
    const base::DPBuffer::WeightsType& weights =
        buffers[time_index].GetWeights();

    for (size_t baseline = 0; baseline < n_baselines_in_buffers; ++baseline) {
      const size_t antenna1 = antennas1[baseline];
      const size_t antenna2 = antennas2[baseline];
      if (antenna1 != antenna2) {
        for (size_t channel_block_index = 0;
             channel_block_index != n_channel_blocks; ++channel_block_index) {
          ChannelBlockData& cb_data = channel_blocks_[channel_block_index];
          size_t& vis_index = visibility_indices[channel_block_index];
          const size_t first_channel = channel_begin[channel_block_index];
          const size_t end_channel = channel_begin[channel_block_index + 1];
          const size_t channel_block_size = end_channel - first_channel;

          for (size_t i = 0; i < channel_block_size; ++i) {
            cb_data.data_[vis_index + i] = ToSpecificMatrix<MatrixType>(
                &data(baseline, first_channel + i, 0));
            cb_data.antenna_indices_[vis_index + i] =
                std::pair<uint32_t, uint32_t>(antenna1, antenna2);
          }
          if (has_weights) {
            for (size_t cb_index = 0; cb_index < channel_block_size;
                 ++cb_index) {
              for (size_t corr = 0; corr != weights.shape()[2]; ++corr)
                cb_data.weights_(vis_index + cb_index, corr) =
                    weights(baseline, first_channel + cb_index, corr);
            }
          }

          for (size_t direction = 0; direction < n_directions; ++direction) {
            const base::DPBuffer::DataType& model_data =
                buffers[time_index].GetData(direction_names[direction]);
            const size_t n_solutions =
                channel_blocks_.front().n_solutions_[direction];
            // Calculate the absolute index as required for solution_map_
            const size_t solution_index = time_index * n_solutions / n_times +
                                          solution_start_indices[direction];

            for (size_t i = 0; i < channel_block_size; ++i) {
              cb_data.model_data_(direction, vis_index + i) =
                  ToSpecificMatrix<MatrixType>(
                      &model_data(baseline, first_channel + i, 0));

              cb_data.solution_map_(direction, vis_index + i) = solution_index;
            }
          }

          vis_index += channel_block_size;
        }
      }
    }
  }

  CountAntennaVisibilities(n_antennas);
}

template <typename MatrixType>
SolveData<MatrixType>::SolveData(const BdaSolverBuffer& buffer,
                                 size_t n_channel_blocks, size_t n_directions,
                                 size_t n_antennas,
                                 const std::vector<int>& antennas1,
                                 const std::vector<int>& antennas2,
                                 bool with_weights)
    : channel_blocks_(n_channel_blocks) {
  // Count nr of visibilities
  std::vector<size_t> counts(n_channel_blocks, 0);
  for (size_t row = 0; row != buffer.GetDataRows().size(); ++row) {
    const BDABuffer::Row& data_row = *buffer.GetDataRows()[row];
    const size_t antenna1 = antennas1[data_row.baseline_nr];
    const size_t antenna2 = antennas2[data_row.baseline_nr];
    if (antenna1 != antenna2) {
      for (size_t channel_block_index = 0;
           channel_block_index != n_channel_blocks; ++channel_block_index) {
        const size_t channel_start =
            channel_block_index * data_row.n_channels / n_channel_blocks;
        const size_t channel_end =
            (channel_block_index + 1) * data_row.n_channels / n_channel_blocks;
        counts[channel_block_index] += channel_end - channel_start;
      }
    }
  }

  // Allocate
  for (size_t cb = 0; cb != n_channel_blocks; ++cb) {
    channel_blocks_[cb].Resize(counts[cb], n_directions);
    if (with_weights) {
      channel_blocks_[cb].ResizeWeights(counts[cb]);
    }
  }

  // Fill
  std::vector<size_t> visibility_indices(n_channel_blocks, 0);
  for (size_t row = 0; row != buffer.GetDataRows().size(); ++row) {
    const BDABuffer::Row& data_row = *buffer.GetDataRows()[row];
    const size_t antenna1 = antennas1[data_row.baseline_nr];
    const size_t antenna2 = antennas2[data_row.baseline_nr];
    if (antenna1 != antenna2) {
      for (size_t channel_block_index = 0;
           channel_block_index != n_channel_blocks; ++channel_block_index) {
        const size_t channel_start =
            channel_block_index * data_row.n_channels / n_channel_blocks;
        const size_t channel_end =
            (channel_block_index + 1) * data_row.n_channels / n_channel_blocks;
        const size_t channel_block_size = channel_end - channel_start;
        size_t& vis_index = visibility_indices[channel_block_index];
        ChannelBlockData& cb_data = channel_blocks_[channel_block_index];

        const std::complex<float>* data_ptr =
            data_row.data + channel_start * data_row.n_correlations;
        for (size_t i = 0; i != channel_block_size; ++i) {
          cb_data.data_[vis_index + i] = ToSpecificMatrix<MatrixType>(
              &data_ptr[i * data_row.n_correlations]);
          cb_data.antenna_indices_[vis_index + i] =
              std::pair<uint32_t, uint32_t>(antenna1, antenna2);
        }
        if (with_weights) {
          for (size_t i = 0; i != channel_block_size; ++i) {
            for (size_t p = 0; p != cb_data.weights_.shape()[2]; ++p) {
              cb_data.weights_(vis_index + i, p) =
                  1.0;  // TODO use real vis weights
            }
          }
        }

        for (size_t dir = 0; dir != n_directions; ++dir) {
          const BDABuffer::Row& model_data_row =
              *buffer.GetModelDataRows(dir)[row];
          const std::complex<float>* model_data_ptr =
              model_data_row.data +
              channel_start * model_data_row.n_correlations;
          for (size_t i = 0; i != channel_block_size; ++i) {
            cb_data.model_data_(dir, vis_index + i) =
                ToSpecificMatrix<MatrixType>(
                    &model_data_ptr[i * model_data_row.n_correlations]);
          }
        }

        vis_index += channel_block_size;
      }
    }
  }

  CountAntennaVisibilities(n_antennas);
  for (ChannelBlockData& cb_data : channel_blocks_)
    cb_data.InitializeSolutionIndices();  // TODO replace!
}

template <typename MatrixType>
void SolveData<MatrixType>::CountAntennaVisibilities(size_t n_antennas) {
  for (ChannelBlockData& cb_data : channel_blocks_) {
    cb_data.antenna_visibility_counts_.assign(n_antennas, 0);
    for (const std::pair<uint32_t, uint32_t>& a : cb_data.antenna_indices_) {
      ++cb_data.antenna_visibility_counts_[a.first];
      ++cb_data.antenna_visibility_counts_[a.second];
    }
  }
}

template <typename MatrixType>
void SolveData<MatrixType>::ChannelBlockData::InitializeSolutionIndices() {
  // This initializes the solution indices for direction-independent intervals
  // TODO support DD intervals
  n_solutions_.assign(NDirections(), 1);
  for (size_t i = 0; i != NDirections(); ++i) {
    xt::view(solution_map_, i, xt::all()).fill(i);
  }
}

template class SolveData<aocommon::MC2x2F>;
template class SolveData<aocommon::MC2x2FDiag>;

}  // namespace ddecal
}  // namespace dp3
