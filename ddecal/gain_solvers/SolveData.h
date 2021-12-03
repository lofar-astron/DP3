// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVE_DATA_H
#define DDECAL_SOLVE_DATA_H

#include <aocommon/matrix2x2.h>

#include <vector>

namespace dp3 {
namespace ddecal {

class BdaSolverBuffer;
class SolverBuffer;

/**
 * Contains exactly the data required for solving: (weighted) data, model_data
 * and the associated antennas for each visibility. In this class, the term
 * visibility refers to a 2x2 matrix containing the 4 polarizations.
 */
class SolveData {
 public:
  class ChannelBlockData {
   public:
    using const_iterator = std::vector<aocommon::MC2x2F>::const_iterator;

    void Resize(size_t n_visibilities, size_t n_directions) {
      data_.resize(n_visibilities);
      model_data_.resize(n_directions);
      for (std::vector<aocommon::MC2x2F>& v : model_data_)
        v.resize(n_visibilities);
      antenna_indices_.resize(n_visibilities);
      n_solutions_.resize(n_directions);
      solution_map_.resize(n_directions);
      for (std::vector<uint32_t>& m : solution_map_) m.resize(n_visibilities);
    }
    size_t NDirections() const { return model_data_.size(); }
    size_t NVisibilities() const { return data_.size(); }
    /***
     * The number of visibilities in which a given antenna participates.
     */
    size_t NAntennaVisibilities(size_t antenna_index) const {
      return antenna_visibility_counts_[antenna_index];
    }
    uint32_t Antenna1Index(size_t visibility_index) const {
      return antenna_indices_[visibility_index].first;
    }
    uint32_t Antenna2Index(size_t visibility_index) const {
      return antenna_indices_[visibility_index].second;
    }
    /**
     * Absolute solution index for a direction and visibility combination.
     * When using direction-dependent intervals, a single direction might
     * have multiple solutions. The solution index is absolute, meaning that
     * direction zero starts counting at zero and indices of subsequent
     * directions start after the previous direction. For every direction, the
     * first solution index is also the lowest value, i.e. SolutionIndex(D, 0)
     * <= SolutionIndex(D, i) for any i,D
     */
    uint32_t SolutionIndex(size_t direction_index,
                           size_t visibility_index) const {
      return solution_map_[direction_index][visibility_index];
    }
    /**
     * Returns a vector that maps visibility indices to solution indices (see
     * @ref SolutionIndex())
     */
    const std::vector<uint32_t>& SolutionMap(size_t direction_index) const {
      return solution_map_[direction_index];
    }

    const aocommon::MC2x2F& Visibility(size_t index) const {
      return data_[index];
    }
    const aocommon::MC2x2F& ModelVisibility(size_t direction,
                                            size_t index) const {
      return model_data_[direction][index];
    }
    const std::vector<aocommon::MC2x2F>& ModelVisibilityVector(
        size_t direction) const {
      return model_data_[direction];
    }

    const_iterator DataBegin() const { return data_.begin(); }
    const_iterator DataEnd() const { return data_.end(); }

    uint32_t NSolutionsForDirection(size_t direction_index) const {
      return n_solutions_[direction_index];
    }

   private:
    friend class SolveData;

    /**
     * Initialize n_solutions_ and solution_map_.
     */
    void InitializeSolutionIndices();

    std::vector<aocommon::MC2x2F> data_;
    // _modelData[D] is the model data for direction D
    std::vector<std::vector<aocommon::MC2x2F>> model_data_;
    // Element i contains the first and second antenna corresponding with
    // _data[i] and _modelData[D][i]
    std::vector<std::pair<uint32_t, uint32_t>> antenna_indices_;
    std::vector<size_t> antenna_visibility_counts_;
    /// number of solutions, indexed by direction
    std::vector<uint32_t> n_solutions_;
    /// solution_map_[D][i] is the solution associated to
    /// direction D, visibility index i.
    std::vector<std::vector<uint32_t>> solution_map_;
  };

  /**
   * Constructor for regular data.
   * @param buffer Buffer with regular data for the current solution interval.
   * @param n_channel_blocks Number of channel blocks / groups.
   * @param n_directions Number of solver directions.
   * @param n_antennas Number of antennas.
   * @param n_solutions_per_direction For each direction, the number of
   * solutions for this direction. The timesteps in the buffer are split evenly
   * over the solutions. This allows direction-dependent solution intervals. If
   * n_solutions_per_direction[i] is larger than the number of available
   * timesteps, it is truncated.
   * @param antennas1 For each baseline, the index of the first antenna.
   * @param antennas2 For each baseline, the index of the second antenna.
   */
  SolveData(const SolverBuffer& buffer, size_t n_channel_blocks,
            size_t n_directions, size_t n_antennas,
            const std::vector<size_t>& n_solutions_per_direction,
            const std::vector<int>& antennas1,
            const std::vector<int>& antennas2);

  /**
   * Constructor for BDA data.
   * @param buffer Buffer with BDA data for the current solution interval.
   * @param n_channel_blocks Number of channel blocks / groups.
   * @param n_directions Number of solver directions.
   * @param n_antennas Number of antennas.
   * @param antennas1 For each baseline, the index of the first antenna.
   * @param antennas2 For each baseline, the index of the second antenna.
   */
  SolveData(const BdaSolverBuffer& buffer, size_t n_channel_blocks,
            size_t n_directions, size_t n_antennas,
            const std::vector<int>& antennas1,
            const std::vector<int>& antennas2);

  size_t NChannelBlocks() const { return channel_blocks_.size(); }

  const ChannelBlockData& ChannelBlock(size_t i) const {
    return channel_blocks_[i];
  }

 private:
  void CountAntennaVisibilities(size_t n_antennas);

  /// The data, indexed by channel block index
  std::vector<ChannelBlockData> channel_blocks_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
