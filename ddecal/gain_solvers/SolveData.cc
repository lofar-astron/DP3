#include "SolveData.h"

#include "BDASolverBuffer.h"

namespace dp3 {
namespace base {

SolveData::SolveData(const BDASolverBuffer& buffer, size_t n_channel_blocks,
                     size_t n_directions, size_t n_antennas,
                     const std::vector<int>& ant1, const std::vector<int>& ant2)
    : channel_blocks_(n_channel_blocks) {
  // Count nr of visibilities
  std::vector<size_t> counts(n_channel_blocks, 0);
  for (size_t row = 0; row != buffer.GetDataRows().size(); ++row) {
    const BDABuffer::Row& data_row = *buffer.GetDataRows()[row];
    const size_t antenna1 = ant1[data_row.baseline_nr];
    const size_t antenna2 = ant2[data_row.baseline_nr];
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
    const size_t antenna1 = ant1[data_row.baseline_nr];
    const size_t antenna2 = ant2[data_row.baseline_nr];
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

        for (size_t i = 0; i != channel_block_size; ++i) {
          cb_data.antenna_indices_[vis_index + i] =
              std::pair<uint32_t, uint32_t>(antenna1, antenna2);
        }

        vis_index += channel_block_size;
      }
    }
  }

  // Count the visibilities per antenna
  for (ChannelBlockData& cb_data : channel_blocks_) {
    cb_data.antenna_visibility_counts_.assign(n_antennas, 0);
    for (const std::pair<uint32_t, uint32_t>& a : cb_data.antenna_indices_) {
      cb_data.antenna_visibility_counts_[a.first]++;
      cb_data.antenna_visibility_counts_[a.second]++;
    }
  }
}

}  // namespace base
}  // namespace dp3
