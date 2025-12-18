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
    size_t solution_interval)
    : n_solutions_per_direction_(n_solutions_per_direction),
      n_sub_solutions_(std::accumulate(n_solutions_per_direction.begin(),
                                       n_solutions_per_direction.end(), 0u)),
      n_directions_(n_solutions_per_direction.size()),
      n_antennas_(n_antennas),
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
    // If e.g. there are [3, 3, 5, 5] solutions per direction, and the solution
    // interval is 15, upsample the solutions to 15 intervals instead of 5.
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
        solutions,
    size_t n_pol) const {
  assert(solutions[0][0].size() == n_antennas_ * n_sub_solutions_ * n_pol);

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
                                                  n_pol);
      for (size_t antenna = 0; antenna != n_antennas_; ++antenna) {
        for (size_t direction = 0; direction != n_directions_; ++direction) {
          for (size_t pol = 0; pol != n_pol; ++pol) {
            // The major time index could be computed in the outer loop,
            // but not yet done for clarity.
            size_t time_index_original;
            size_t offset_original;
            std::tie(time_index_original, offset_original) =
                MapResampledToOriginal(time_index, antenna, direction, pol,
                                       n_pol);
            const size_t offset_upsampled =
                (antenna * n_directions_ + direction) * n_pol + pol;
            upsampled[time_index][channel_block][offset_upsampled] =
                solutions[time_index_original][channel_block][offset_original];
          }
        }
      }
    }
  }
  return upsampled;
}

std::vector<std::vector<std::vector<ConstraintResult>>>
SolutionResampler::Upsample(
    const std::vector<std::vector<std::vector<ConstraintResult>>>& solutions)
    const {
  const size_t n_intervals = solutions.size();
  const size_t n_constraints = n_intervals == 0 ? 0 : solutions.front().size();
  const size_t n_constraint_sub_types =
      n_constraints == 0 ? 0 : solutions.front().front().size();
  std::vector<std::vector<std::vector<ConstraintResult>>> upsampled;
  upsampled.resize(n_intervals * n_substeps_);

  // Fill the upsampled array with meta data and allocate data and weight arrays
  for (size_t interval_index = 0; interval_index != n_intervals;
       ++interval_index) {
    const size_t new_interval_offset = interval_index * n_substeps_;
    for (size_t sub_solution = 0; sub_solution != n_substeps_; ++sub_solution) {
      size_t new_interval_index = new_interval_offset + sub_solution;
      upsampled[new_interval_index].resize(n_constraints);
      for (size_t constraint_index = 0; constraint_index != n_constraints;
           ++constraint_index) {
        upsampled[new_interval_index][constraint_index].resize(
            n_constraint_sub_types);
        for (size_t sub_type_index = 0;
             sub_type_index != n_constraint_sub_types; ++sub_type_index) {
          const ConstraintResult& input_result =
              solutions[interval_index][constraint_index][sub_type_index];
          const size_t n_input_values = input_result.vals.size();
          assert(input_result.weights.size() == n_input_values);
          assert(input_result.dims[1] == n_sub_solutions_);
          const size_t n_new_values =
              n_input_values * n_directions_ / n_sub_solutions_;

          ConstraintResult& upsampled_result =
              upsampled[new_interval_index][constraint_index][sub_type_index];
          upsampled_result.vals.resize(n_new_values);
          upsampled_result.weights.resize(n_new_values);
          upsampled_result.axes = input_result.axes;
          upsampled_result.dims = input_result.dims;
          upsampled_result.dims[1] = n_directions_;
          upsampled_result.name = input_result.name;
        }
      }
    }
  }

  // Do the upsampling
  for (size_t constraint_index = 0; constraint_index != n_constraints;
       ++constraint_index) {
    for (size_t sub_type_index = 0; sub_type_index != n_constraint_sub_types;
         ++sub_type_index) {
      for (size_t interval_index = 0; interval_index != n_intervals;
           ++interval_index) {
        const ConstraintResult& input_result =
            solutions[interval_index][constraint_index][sub_type_index];
        const size_t n_input_values = input_result.vals.size();
        const size_t new_interval_offset = interval_index * n_substeps_;

        size_t old_sub_offset = 0;
        for (size_t direction = 0; direction != n_directions_; ++direction) {
          const size_t n_sub_solutions = n_solutions_per_direction_[direction];
          for (size_t sub_solution = 0; sub_solution != n_substeps_;
               ++sub_solution) {
            const size_t old_sub_index =
                old_sub_offset + sub_solution * n_sub_solutions / n_substeps_;
            const size_t new_interval_index =
                new_interval_offset + sub_solution;
            ConstraintResult& destination =
                upsampled[new_interval_index][constraint_index][sub_type_index];

            // Copy the data for this substep and direction. In the constraints,
            // the subsolution axis is always dimension 1. Dimension 0 is the
            // antenna axis. There may be more axes for e.g. frequency (often
            // dimension 2). Antenna is increasing slowest.
            const size_t dim0 = input_result.dims[0];
            const size_t residual_dims =
                n_input_values / (dim0 * n_sub_solutions_);
            for (size_t i = 0; i != dim0; ++i) {
              const size_t in_index =
                  (i * n_sub_solutions_ + old_sub_index) * residual_dims;
              const size_t out_index =
                  (i * n_directions_ + direction) * residual_dims;
              std::copy_n(input_result.vals.begin() + in_index, residual_dims,
                          destination.vals.begin() + out_index);
              std::copy_n(input_result.weights.begin() + in_index,
                          residual_dims,
                          destination.weights.begin() + out_index);
            }
          }
          old_sub_offset += n_sub_solutions;
        }
      }  // loop over interval_index
    }    // loop over sub_type_index
  }      // loop over constraint_index
  return upsampled;
}

std::pair<size_t, size_t> SolutionResampler::MapResampledToOriginal(
    size_t time_index, size_t antenna_index, size_t direction_index,
    size_t pol_index, size_t n_pol) const {
  const size_t major_time_index = time_index / n_substeps_;
  size_t substep_index = time_index % n_substeps_;

  const size_t antenna_offset = antenna_index * n_sub_solutions_ * n_pol;
  const size_t direction_offset =
      (std::accumulate(n_solutions_per_direction_.begin(),
                       n_solutions_per_direction_.begin() + direction_index,
                       0u) +
       substep_index /
           (n_substeps_ / n_solutions_per_direction_[direction_index])) *
      n_pol;

  const size_t inner_offset = antenna_offset + direction_offset + pol_index;
  return std::make_pair(major_time_index, inner_offset);
}

}  // namespace ddecal
}  // namespace dp3
