// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "KernelSmoother.h"
#include "SmoothnessConstraint.h"
#include <aocommon/parallelfor.h>
#include <boost/make_unique.hpp>

namespace dp3 {
namespace ddecal {

SmoothnessConstraint::SmoothnessConstraint(double bandwidth_hz,
                                           double bandwidth_ref_frequency_hz)
    : kernel_type_(Smoother::GaussianKernel),
      bandwidth_(bandwidth_hz),
      bandwidth_ref_frequency_(bandwidth_ref_frequency_hz) {}

void SmoothnessConstraint::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);
  frequencies_ = frequencies;
  for (int32_t v : solutions_per_direction) {
    if (v != 1)
      throw std::runtime_error(
          "The smoothness constraint does not yet support direction-dependent "
          "intervals");
  }
}

void SmoothnessConstraint::SetDistanceFactors(
    std::vector<double>&& antenna_distance_factors) {
  antenna_distance_factors_ = std::move(antenna_distance_factors);
  if (!loop_) {
    loop_ = boost::make_unique<aocommon::ParallelFor<size_t>>(NThreads());
  }
  for (size_t i = 0; i != NThreads(); ++i)
    fit_data_.emplace_back(frequencies_, kernel_type_, bandwidth_,
                           bandwidth_ref_frequency_);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  const size_t n_pol = solutions.front().size() / (NAntennas() * NDirections());

  loop_->Run(
      0, NAntennas() * NDirections(), [&](size_t ant_dir_index, size_t thread) {
        size_t ant_index = ant_dir_index / NDirections();
        for (size_t pol = 0; pol != n_pol; ++pol) {
          size_t solution_index = ant_dir_index * n_pol + pol;
          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            // Flag channels where calibration yielded inf or nan
            if (isfinite(solutions[ch][solution_index])) {
              fit_data_[thread].data[ch] = solutions[ch][solution_index];
              fit_data_[thread].weight[ch] =
                  weights_[ant_index * NChannelBlocks() + ch];
            } else {
              fit_data_[thread].data[ch] = 0.0;
              fit_data_[thread].weight[ch] = 0.0;
            }
          }

          fit_data_[thread].smoother.Smooth(
              fit_data_[thread].data.data(), fit_data_[thread].weight.data(),
              antenna_distance_factors_[ant_index]);

          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            solutions[ch][solution_index] = fit_data_[thread].data[ch];
          }
        }
      });

  return std::vector<Constraint::Result>();
}

}  // namespace ddecal
}  // namespace dp3
