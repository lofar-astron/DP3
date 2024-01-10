// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVE_DATA_H
#define DDECAL_SOLVE_DATA_H

#include <cstddef>
#include <vector>

#include <aocommon/matrix2x2.h>
#include <xtensor/xtensor.hpp>

#include <dp3/base/DPBuffer.h>

namespace dp3 {
namespace ddecal {

class BdaSolverBuffer;

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
      model_data_.resize({n_directions, n_visibilities});
      antenna_indices_.resize(n_visibilities);
      n_solutions_.resize(n_directions);
      solution_map_.resize({n_directions, n_visibilities});
    }
    void ResizeWeights(size_t n_visibilities) {
      constexpr size_t kNCorrelations = 4;
      weights_.resize({n_visibilities, kNCorrelations});
    }
    size_t NDirections() const { return model_data_.shape(0); }
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
     * <= SolutionIndex(D, i) for any i,D.
     */
    uint32_t SolutionIndex(size_t direction_index,
                           size_t visibility_index) const {
      return solution_map_(direction_index, visibility_index);
    }

    const uint32_t* SolutionMapData() const { return solution_map_.data(); }

    const float& Weight(size_t index) const { return weights_(index, 0); }
    const aocommon::MC2x2F& Visibility(size_t index) const {
      return data_[index];
    }
    const aocommon::MC2x2F& ModelVisibility(size_t direction,
                                            size_t index) const {
      return model_data_(direction, index);
    }

    const_iterator DataBegin() const { return data_.begin(); }
    const_iterator DataEnd() const { return data_.end(); }

    /**
     * @return The number of solutions for the given direction.
     */
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
    // weights_(i, pol) contains the weight for data_[i][pol]. The vector will
    // be left empty when the algorithm does not need the weights.
    xt::xtensor<float, 2> weights_;
    // model_data_(d, i) is the model data for direction d, element i
    xt::xtensor<aocommon::MC2x2F, 2> model_data_;
    // Element i contains the first and second antenna corresponding with
    // data_[i] and model_data_(d, i)
    std::vector<std::pair<uint32_t, uint32_t>> antenna_indices_;
    std::vector<size_t> antenna_visibility_counts_;
    /// number of solutions, indexed by direction
    std::vector<uint32_t> n_solutions_;
    /// solution_map_(D,i) is the solution associated to
    /// direction D, visibility index i.
    xt::xtensor<uint32_t, 2> solution_map_;
  };

  /**
   * Constructor for regular data.
   * @param buffers Weighted data and weighted model data for all time steps in
   * the current solution interval.
   * @param directions_names Names of the model data in 'buffers'.
   * @param n_channel_blocks Number of channel blocks / groups.
   * @param n_antennas Number of antennas.
   * @param n_solutions_per_direction For each direction, the number of
   * solutions for this direction. The timesteps in the buffer are split evenly
   * over the solutions. This allows direction-dependent solution intervals. If
   * n_solutions_per_direction[i] is larger than the number of available
   * timesteps, it is truncated.
   * @param antennas1 For each baseline, the index of the first antenna.
   * @param antennas2 For each baseline, the index of the second antenna.
   */
  SolveData(const std::vector<base::DPBuffer>& buffers,
            const std::vector<std::string>& direction_names,
            size_t n_channel_blocks, size_t n_antennas,
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
            const std::vector<int>& antennas2, bool with_weights);

  size_t NChannelBlocks() const { return channel_blocks_.size(); }

  const ChannelBlockData& ChannelBlock(size_t i) const {
    return channel_blocks_[i];
  }

 private:
  void CountAntennaVisibilities(size_t n_antennas);

  /// The data, indexed by channel block index
  std::vector<ChannelBlockData> channel_blocks_;
};

template <bool Add>
void DiagonalAddOrSubtractDirection(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    size_t n_solutions, const std::vector<std::complex<double>>& solutions) {
  using DComplex = std::complex<double>;
  using Complex = std::complex<float>;
  using aocommon::MC2x2F;
  const size_t n_visibilities = cb_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
    const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
    const uint32_t solution_index = cb_data.SolutionIndex(direction, vis_index);
    const DComplex* solution_1 =
        &solutions[(antenna_1 * n_solutions + solution_index) * 2];
    const DComplex* solution_2 =
        &solutions[(antenna_2 * n_solutions + solution_index) * 2];
    const Complex solution_1_0(solution_1[0]);
    const Complex solution_1_1(solution_1[1]);
    const Complex solution_2_0_conj(std::conj(solution_2[0]));
    const Complex solution_2_1_conj(std::conj(solution_2[1]));

    MC2x2F& data = v_residual[vis_index];
    const MC2x2F& model = cb_data.ModelVisibility(direction, vis_index);
    const MC2x2F contribution(solution_1_0 * model[0] * solution_2_0_conj,
                              solution_1_0 * model[1] * solution_2_1_conj,
                              solution_1_1 * model[2] * solution_2_0_conj,
                              solution_1_1 * model[3] * solution_2_1_conj);
    if constexpr (Add)
      data += contribution;
    else
      data -= contribution;
  }
}

}  // namespace ddecal
}  // namespace dp3

#endif
