// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDE_SOLVE_DATA_H
#define DDE_SOLVE_DATA_H

#include <aocommon/matrix2x2.h>

#include <vector>

namespace dp3 {
namespace ddecal {

class BDASolverBuffer;

/**
 * Contains exactly the data required for solving: (weighted) data, model_data
 * and the associated antennas for each visibility. In this class, the term
 * visibility refers to a 2x2 matrix containing the 4 polarizations.
 */
class SolveData {
 public:
  class ChannelBlockData {
   public:
    void Resize(size_t n_visibilities, size_t n_directions) {
      data_.resize(n_visibilities);
      model_data_.resize(n_directions);
      for (std::vector<aocommon::MC2x2F>& v : model_data_)
        v.resize(n_visibilities);
      antenna_indices_.resize(n_visibilities);
    }
    size_t NDirections() const { return model_data_.size(); }
    size_t NVisibilities() const { return data_.size(); }
    /***
     * The number of visibilities in which a given antenna participates.
     */
    size_t NAntennaVisibilities(size_t antenna_index) const {
      return antenna_visibility_counts_[antenna_index];
    }
    size_t Antenna1Index(size_t visibility_index) const {
      return antenna_indices_[visibility_index].first;
    }
    size_t Antenna2Index(size_t visibility_index) const {
      return antenna_indices_[visibility_index].second;
    }

    const aocommon::MC2x2F& Visibility(size_t index) const {
      return data_[index];
    }
    const aocommon::MC2x2F& ModelVisibility(size_t dir, size_t index) const {
      return model_data_[dir][index];
    }

   private:
    friend class SolveData;
    std::vector<aocommon::MC2x2F> data_;
    // _modelData[D] is the model data for direction D
    std::vector<std::vector<aocommon::MC2x2F>> model_data_;
    // Element i contains the first and second antenna corresponding with
    // _data[i] and _modelData[D][i]
    std::vector<std::pair<uint32_t, uint32_t>> antenna_indices_;
    std::vector<size_t> antenna_visibility_counts_;
  };

  SolveData(const BDASolverBuffer& buffer, size_t n_channel_blocks,
            size_t n_directions, size_t n_antennas,
            const std::vector<int>& ant1, const std::vector<int>& ant2);

  size_t NChannelBlocks() const { return channel_blocks_.size(); }

  ChannelBlockData& ChannelBlock(size_t i) { return channel_blocks_[i]; }
  const ChannelBlockData& ChannelBlock(size_t i) const {
    return channel_blocks_[i];
  }

 private:
  std::vector<ChannelBlockData> channel_blocks_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
