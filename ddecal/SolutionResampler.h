// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLUTIONRESAMPLER_H_
#define DP3_DDECAL_SOLUTIONRESAMPLER_H_

#include <complex>
#include <vector>

namespace dp3 {
namespace ddecal {

/**
 * @brief Class for resampling the solution into a square multi-dimensional
 * array.
 */
class SolutionResampler {
 public:
  SolutionResampler(const std::vector<size_t>& n_solutions_per_direction,
                    size_t n_antennas, size_t n_pol, size_t solution_interval);
  ~SolutionResampler() = default;

  /**
   * @brief Compute the upsampled solutions for dd interval solutions.
   *
   * Dimension of returned vector is [n_times * n_substeps, n_channel_blocks,
   * n_antenna * n_directions * n_pol]
   *
   * @param solutions Original solutions
   * @return std::vector<std::vector<std::vector<std::complex<double>>>>
   * Upsampled solutions
   */
  std::vector<std::vector<std::vector<std::complex<double>>>> Upsample(
      const std::vector<std::vector<std::vector<std::complex<double>>>>&
          solutions) const;

  /**
   * @brief Get the number of substeps per solution interval.
   *
   *  NOTE: the current implementation for computing the number of substeps
   * assumes that the solution interval is an integer multiple of the values in
   * the number of solutions per direction.
   */
  size_t GetNrSubSteps() const { return n_substeps_; }

  /**
   * @brief Map index in upsampled solutions vector to the matching
   * index in the original solution vector.
   */
  std::pair<size_t, size_t> MapResampledToOriginal(size_t time_index,
                                                   size_t antenna_index,
                                                   size_t direction_index,
                                                   size_t pol_index) const;

 private:
  const std::vector<size_t> n_solutions_per_direction_;
  const size_t n_solutions_;
  const size_t n_directions_;

  const size_t n_antennas_;
  const size_t n_pol_;
  const size_t solution_interval_;
  size_t n_substeps_;
};

}  // namespace ddecal
}  // namespace dp3

#endif