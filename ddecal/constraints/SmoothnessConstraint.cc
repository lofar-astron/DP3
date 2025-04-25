// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SmoothnessConstraint.h"

#include <aocommon/staticfor.h>

namespace dp3 {
namespace ddecal {

SmoothnessConstraint::SmoothnessConstraint(double bandwidth_hz,
                                           double bandwidth_ref_frequency_hz,
                                           double spectral_exponent,
                                           bool kernel_truncation)
    : kernel_type_(Smoother::GaussianKernel),
      bandwidth_(bandwidth_hz),
      bandwidth_ref_frequency_(bandwidth_ref_frequency_hz),
      spectral_exponent_(spectral_exponent),
      kernel_truncation_(kernel_truncation) {}

void SmoothnessConstraint::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);
  frequencies_ = frequencies;
  for (size_t i = 0; i != aocommon::ThreadPool::GetInstance().NThreads(); ++i) {
    fit_data_.emplace_back(frequencies_, kernel_type_, bandwidth_,
                           bandwidth_ref_frequency_, spectral_exponent_,
                           kernel_truncation_);
  }
}

void SmoothnessConstraint::SetDdSmoothingFactors(
    std::vector<double> dd_smoothing_factors) {
  if (!dd_smoothing_factors.empty() &&
      dd_smoothing_factors.size() != NSubSolutions()) {
    throw std::runtime_error(
        "Invalid setting for the dd factors option of smoothness constraint: "
        "the specified factors is a vector of size " +
        std::to_string(dd_smoothing_factors.size()) +
        ", whereas the number of solutions is " +
        std::to_string(NSubSolutions()));
  }
  dd_smoothing_factors_ = std::move(dd_smoothing_factors);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    SolutionSpan& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  assert(NChannelBlocks() == solutions.shape(0));
  assert(NAntennas() == solutions.shape(1));
  assert(NSubSolutions() == solutions.shape(2));
  assert(dd_smoothing_factors_.empty() ||
         NSubSolutions() == dd_smoothing_factors_.size());

  const size_t n_polarizations = solutions.shape(3);
  const size_t n_smoothed = NAntennas() * NSubSolutions() * n_polarizations;
  auto solutions_view =
      xt::reshape_view(solutions, {NChannelBlocks(), n_smoothed});

  aocommon::StaticFor<size_t> loop;
  loop.Run(
      0, n_smoothed, [&](size_t begin_index, size_t end_index, size_t thread) {
        for (size_t smoothing_index = begin_index; smoothing_index < end_index;
             ++smoothing_index) {
          const size_t sol_index =
              (smoothing_index / n_polarizations) % NSubSolutions();
          const size_t ant_index =
              smoothing_index / (NSubSolutions() * n_polarizations);
          const double* weights = sub_solution_weights_.empty()
                                      ? weights_.data()
                                      : sub_solution_weights_[sol_index].data();
          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            // Flag channels where calibration yielded inf or nan
            if (isfinite(solutions_view(ch, smoothing_index))) {
              fit_data_[thread].data[ch] = solutions_view(ch, smoothing_index);
              fit_data_[thread].weight[ch] =
                  weights[ant_index * NChannelBlocks() + ch];
            } else {
              fit_data_[thread].data[ch] = 0.0;
              fit_data_[thread].weight[ch] = 0.0;
            }
          }

          const double dd_factor = dd_smoothing_factors_.empty()
                                       ? 1.0
                                       : 1.0 / dd_smoothing_factors_[sol_index];

          fit_data_[thread].smoother.Smooth(
              fit_data_[thread].data.data(), fit_data_[thread].weight.data(),
              antenna_factors_[ant_index] * dd_factor);

          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            solutions_view(ch, smoothing_index) = fit_data_[thread].data[ch];
          }
        }
      });

  return std::vector<Constraint::Result>();
}

}  // namespace ddecal
}  // namespace dp3
