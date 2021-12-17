// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionResampler.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <tuple>

namespace dp3 {
namespace ddecal {

SolutionResampler::SolutionResampler(
    const std::vector<size_t>& n_solutions_per_direction, size_t n_antennas,
    size_t n_pol, size_t solution_interval)
    : n_solutions_per_direction_(n_solutions_per_direction),
      n_solutions_(std::accumulate(n_solutions_per_direction.begin(),
                                   n_solutions_per_direction.end(), 0u)),
      n_directions_(n_solutions_per_direction.size()),
      n_antennas_(n_antennas),
      n_pol_(n_pol),
      solution_interval_(solution_interval) {
  // Compute the number of substeps per solution interval.
  // Note that the computed nr_substeps might be a conservative approximation
  // for the least common multiple.
  const size_t max_n_solutions = *std::max_element(
      n_solutions_per_direction_.begin(), n_solutions_per_direction_.end());
  bool use_solution_interval = false;
  for (const size_t val : n_solutions_per_direction_) {
    // Basic assumption: values in n_solutions_per_direction are integer
    // divisors of solution_interval
    assert(solution_interval_ % val == 0);
    if (max_n_solutions % val != 0) {
      use_solution_interval = true;
      break;
    }
  }
  n_substeps_ = use_solution_interval ? solution_interval : max_n_solutions;
}

std::vector<std::vector<std::vector<std::complex<double>>>>
SolutionResampler::Upsample(
    const std::vector<std::vector<std::vector<std::complex<double>>>>&
        solutions) const {
  assert(solutions[0][0].size() == n_antennas_ * n_solutions_ * n_pol_);

  // This assumes the number of channel blocks is identical for each time
  const size_t n_channel_blocks = solutions[0].size();

  // Initialize outer dimension of upsampled solution
  std::vector<std::vector<std::vector<std::complex<double>>>> upsampled;
  upsampled.resize(solutions.size() * n_substeps_);

  for (size_t time_index = 0; time_index != upsampled.size(); ++time_index) {
    upsampled[time_index].resize(n_channel_blocks);
    for (size_t channel_block = 0; channel_block != n_channel_blocks;
         ++channel_block) {
      upsampled[time_index][channel_block].resize(n_antennas_ * n_directions_ *
                                                  n_pol_);
      for (size_t antenna = 0; antenna != n_antennas_; ++antenna) {
        for (size_t direction = 0; direction != n_directions_; ++direction) {
          for (size_t pol = 0; pol != n_pol_; ++pol) {
            // The major time index could be computed in the outer loop,
            // but not yet done for clarity.
            size_t time_index_original;
            size_t offset_original;
            std::tie(time_index_original, offset_original) =
                MapResampledToOriginal(time_index, antenna, direction, pol);
            const size_t offset_upsampled =
                (antenna * n_directions_ + direction) * n_pol_ + pol;
            upsampled[time_index][channel_block][offset_upsampled] =
                solutions[time_index_original][channel_block][offset_original];
          }
        }
      }
    }
  }
  return upsampled;
}

std::pair<size_t, size_t> SolutionResampler::MapResampledToOriginal(
    size_t time_index, size_t antenna_index, size_t direction_index,
    size_t pol_index) const {
  const size_t major_time_index = time_index / n_substeps_;
  size_t substep_index = time_index % n_substeps_;

  const size_t antenna_offset = antenna_index * n_solutions_ * n_pol_;
  const size_t direction_offset =
      (std::accumulate(n_solutions_per_direction_.begin(),
                       n_solutions_per_direction_.begin() + direction_index,
                       0u) +
       substep_index /
           (n_substeps_ / n_solutions_per_direction_[direction_index])) *
      n_pol_;

  const size_t inner_offset = antenna_offset + direction_offset + pol_index;
  return std::make_pair(major_time_index, inner_offset);
}

}  // namespace ddecal
}  // namespace dp3