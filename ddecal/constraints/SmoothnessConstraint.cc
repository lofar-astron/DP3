// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "KernelSmoother.h"
#include "SmoothnessConstraint.h"
#include <aocommon/dynamicfor.h>

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
}

void SmoothnessConstraint::SetDistanceFactors(
    std::vector<double>&& antenna_distance_factors) {
  antenna_distance_factors_ = std::move(antenna_distance_factors);
  if (!loop_) {
    loop_ = std::make_unique<aocommon::DynamicFor<size_t>>(NThreads());
  }
  for (size_t i = 0; i != NThreads(); ++i)
    fit_data_.emplace_back(frequencies_, kernel_type_, bandwidth_,
                           bandwidth_ref_frequency_);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    SolutionSpan& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  assert(NChannelBlocks() == solutions.shape(0));
  assert(NAntennas() == solutions.shape(1));
  assert(NSolutions() == solutions.shape(2));
  const size_t n_polarizations = solutions.shape(3);
  const size_t n_smoothed = NAntennas() * NSolutions() * n_polarizations;
  auto solutions_view =
      xt::reshape_view(solutions, {NChannelBlocks(), n_smoothed});

  loop_->Run(0, n_smoothed, [&](size_t smoothing_index, size_t thread) {
    size_t ant_index = smoothing_index / (NSolutions() * n_polarizations);
    for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
      // Flag channels where calibration yielded inf or nan
      if (isfinite(solutions_view(ch, smoothing_index))) {
        fit_data_[thread].data[ch] = solutions_view(ch, smoothing_index);
        fit_data_[thread].weight[ch] =
            weights_[ant_index * NChannelBlocks() + ch];
      } else {
        fit_data_[thread].data[ch] = 0.0;
        fit_data_[thread].weight[ch] = 0.0;
      }
    }

    fit_data_[thread].smoother.Smooth(fit_data_[thread].data.data(),
                                      fit_data_[thread].weight.data(),
                                      antenna_distance_factors_[ant_index]);

    for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
      solutions_view(ch, smoothing_index) = fit_data_[thread].data[ch];
    }
  });

  return std::vector<Constraint::Result>();
}

}  // namespace ddecal
}  // namespace dp3
