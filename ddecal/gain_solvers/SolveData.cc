// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolveData.h"

#include "BdaSolverBuffer.h"
#include "SolverBuffer.h"

#include <cassert>

using dp3::base::BDABuffer;

namespace dp3 {
namespace ddecal {

SolveData::SolveData(const SolverBuffer& buffer, size_t n_channel_blocks,
                     size_t n_directions, size_t n_antennas,
                     const std::vector<int>& antennas1,
                     const std::vector<int>& antennas2)
    : channel_blocks_(n_channel_blocks) {
  std::vector<size_t> channel_begin(n_channel_blocks + 1, 0);

  // Count nr of baselines with different antennas.
  size_t n_baselines = 0;
  for (size_t baseline = 0; baseline < buffer.NBaselines(); ++baseline) {
    assert(size_t(antennas1[baseline]) < n_antennas &&
           size_t(antennas2[baseline]) < n_antennas);
    if (antennas1[baseline] != antennas2[baseline]) ++n_baselines;
  }

  // Count nr of visibilities per channel block and allocate memory.
  for (size_t channel_block_index = 0; channel_block_index != n_channel_blocks;
       ++channel_block_index) {
    ChannelBlockData& cb_data = channel_blocks_[channel_block_index];

    channel_begin[channel_block_index + 1] =
        (channel_block_index + 1) * buffer.NChannels() / n_channel_blocks;
    const size_t channel_block_size = channel_begin[channel_block_index + 1] -
                                      channel_begin[channel_block_index];

    cb_data.Resize(buffer.NTimes() * n_baselines * channel_block_size,
                   n_directions);
  }

  // Fill all channel blocks with data.
  std::vector<size_t> visibility_indices(n_channel_blocks, 0);
  for (size_t time_index = 0; time_index < buffer.NTimes(); ++time_index) {
    for (size_t baseline = 0; baseline < buffer.NBaselines(); ++baseline) {
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
            cb_data.data_[vis_index + i] = aocommon::MC2x2F(
                buffer.DataPointer(time_index, baseline, first_channel + i));
            cb_data.antenna_indices_[vis_index + i] =
                std::pair<uint32_t, uint32_t>(antenna1, antenna2);
          }

          for (size_t direction = 0; direction < n_directions; ++direction) {
            for (size_t i = 0; i < channel_block_size; ++i) {
              cb_data.model_data_[direction][vis_index + i] =
                  aocommon::MC2x2F(buffer.ModelDataPointer(
                      time_index, direction, baseline, first_channel + i));
            }
          }

          vis_index += channel_block_size;
        }
      }
    }
  }

  CountAntennaVisibilities(n_antennas);
}

SolveData::SolveData(const BdaSolverBuffer& buffer, size_t n_channel_blocks,
                     size_t n_directions, size_t n_antennas,
                     const std::vector<int>& antennas1,
                     const std::vector<int>& antennas2)
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
          cb_data.data_[vis_index + i] =
              aocommon::MC2x2F(&data_ptr[i * data_row.n_correlations]);
          cb_data.antenna_indices_[vis_index + i] =
              std::pair<uint32_t, uint32_t>(antenna1, antenna2);
        }

        for (size_t dir = 0; dir != n_directions; ++dir) {
          const BDABuffer::Row& model_data_row =
              *buffer.GetModelDataRows(dir)[row];
          const std::complex<float>* model_data_ptr =
              model_data_row.data +
              channel_start * model_data_row.n_correlations;
          for (size_t i = 0; i != channel_block_size; ++i) {
            cb_data.model_data_[dir][vis_index + i] = aocommon::MC2x2F(
                &model_data_ptr[i * model_data_row.n_correlations]);
          }
        }

        vis_index += channel_block_size;
      }
    }
  }

  CountAntennaVisibilities(n_antennas);
}

void SolveData::CountAntennaVisibilities(size_t n_antennas) {
  for (ChannelBlockData& cb_data : channel_blocks_) {
    cb_data.antenna_visibility_counts_.assign(n_antennas, 0);
    for (const std::pair<uint32_t, uint32_t>& a : cb_data.antenna_indices_) {
      ++cb_data.antenna_visibility_counts_[a.first];
      ++cb_data.antenna_visibility_counts_[a.second];
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
