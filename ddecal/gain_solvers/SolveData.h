// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVE_DATA_H
#define DDECAL_SOLVE_DATA_H

#include <cstddef>
#include <numeric>
#include <vector>

#include <aocommon/matrix2x2.h>
#include <aocommon/staticfor.h>
#include <xtensor/xtensor.hpp>

#include "base/DPBuffer.h"

namespace dp3 {
namespace ddecal {

class BdaSolverBuffer;

/**
 * Contains exactly the data required for solving: (weighted) data, model_data
 * and the associated antennas for each visibility. In this class, the term
 * visibility refers to either a 2x2 diagonal or a full 2x2 matrix, containing 2
 * or 4 polarizations respectively.
 */
template <typename MatrixType = aocommon::MC2x2F>
class SolveData {
 public:
  class ChannelBlockData {
   public:
    using const_iterator = typename std::vector<MatrixType>::const_iterator;

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
    /**
     * The number of visibilities in which a given antenna participates.
     */
    size_t NAntennaVisibilities(size_t antenna_index) const {
      return antenna_visibility_counts_[antenna_index];
    }
    /**
     * The number of visibilities in which a given solution participates.
     * This is thus the number of visibilities for one direction and
     * (within that direction) one solution interval, which in the case
     * of regular data means it is the number of channels times the
     * number of baselines times the number of timesteps within the
     * given solution interval. The solution_index should refer to a
     * solution interval that is inside the given direction.
     */
    size_t NSolutionVisibilities(size_t direction_index,
                                 uint32_t solution_index) const;

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
    const MatrixType& Visibility(size_t index) const { return data_[index]; }
    const MatrixType& ModelVisibility(size_t direction, size_t index) const {
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

    /**
     * @return Total number of solutions. This is the sum of the solutions
     * over all directions.
     */
    uint32_t NSubSolutions() const {
      return std::accumulate(n_solutions_.begin(), n_solutions_.end(),
                             uint32_t(0));
    }

   private:
    friend class SolveData<MatrixType>;

    std::vector<MatrixType> data_;
    // weights_(i, pol) contains the weight for data_[i][pol]. The vector will
    // be left empty when the algorithm does not need the weights.
    xt::xtensor<float, 2> weights_;
    // model_data_(d, i) is the model data for direction d, element i
    xt::xtensor<MatrixType, 2> model_data_;
    // Element i contains the first and second antenna corresponding with
    // data_[i] and model_data_(d, i). Using uint32_t instead of size_t reduces
    // memory usage and improves cache/memory performance (same holds for
    // n_solutions_ and solution_map_).
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
            size_t n_antennas,
            const std::vector<size_t>& n_solutions_per_direction,
            const std::vector<int>& antennas1,
            const std::vector<int>& antennas2, bool with_weights);

  size_t NChannelBlocks() const { return channel_blocks_.size(); }

  const ChannelBlockData& ChannelBlock(size_t i) const {
    return channel_blocks_[i];
  }

  /**
   * Get solution weights, which are the direction-dependent weights.
   * @returns n_solution vectors, each of which is an
   * n_antennas * n_channel_blocks vector, where the channel
   * index varies fastest.
   * The total weight is the sum of the absolute value of all visibilities,
   * i.e. the L_1 norm.
   */
  std::vector<std::vector<double>> GetSolutionWeights() const;

 private:
  void CountAntennaVisibilities();

  /// The data, indexed by channel block index
  std::vector<ChannelBlockData> channel_blocks_;
  size_t n_antennas_;
};

/// Stores all 4 polarizations of the data.
using FullSolveData = SolveData<aocommon::MC2x2F>;
/// Stores only the 2 diagonal values of the data (e.g. XX/YY). Because the word
/// "diagonal" is extensively used for the solution type, the term "Duo" is used
/// for this.
using DuoSolveData = SolveData<aocommon::MC2x2FDiag>;
using UniSolveData = SolveData<std::complex<float>>;

template <bool Add, typename MatrixType>
void DiagonalAddOrSubtractDirection(
    const typename SolveData<MatrixType>::ChannelBlockData& cb_data,
    std::vector<MatrixType>& v_residual, size_t direction, size_t n_solutions,
    const std::vector<std::complex<double>>& solutions, size_t max_threads) {
  using DComplex = std::complex<double>;
  using Complex = std::complex<float>;
  const size_t n_visibilities = cb_data.NVisibilities();
  aocommon::RunConstrainedStaticFor<size_t>(
      0, n_visibilities, max_threads,
      [&](size_t start_vis_index, size_t end_vis_index) {
        for (size_t vis_index = start_vis_index; vis_index != end_vis_index;
             ++vis_index) {
          const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
          const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
          const uint32_t solution_index =
              cb_data.SolutionIndex(direction, vis_index);
          const DComplex* solution_1 =
              &solutions[(antenna_1 * n_solutions + solution_index) * 2];
          const DComplex* solution_2 =
              &solutions[(antenna_2 * n_solutions + solution_index) * 2];

          MatrixType& data = v_residual[vis_index];
          const MatrixType& model =
              cb_data.ModelVisibility(direction, vis_index);
          // Convert first to single precision to make calculation easier
          const aocommon::MC2x2FDiag solution_a(
              static_cast<Complex>(solution_1[0]),
              static_cast<Complex>(solution_1[1]));
          const aocommon::MC2x2FDiag solution_b(
              static_cast<Complex>(solution_2[0]),
              static_cast<Complex>(solution_2[1]));
          const MatrixType contribution(solution_a * model *
                                        HermTranspose(solution_b));
          if constexpr (Add)
            data += contribution;
          else
            data -= contribution;
        }
      });
}

extern template class SolveData<std::complex<float>>;
extern template class SolveData<aocommon::MC2x2F>;
extern template class SolveData<aocommon::MC2x2FDiag>;

}  // namespace ddecal
}  // namespace dp3

#endif
