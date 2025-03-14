// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_
#define DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_

#include <vector>

#include "Constraint.h"
#include "KernelSmoother.h"

namespace dp3 {
namespace ddecal {

class SmoothnessConstraint final : public Constraint {
 public:
  /**
   * @param bandwidth_hz Size of the kernel (smoothing strength)
   * @param bandwidth_ref_frequency_hz Reference frequency for the kernel size,
   * may be zero to have a constant kernel size over frequency.
   * @sa KernelSmoother documentation.
   */
  SmoothnessConstraint(double bandwidth_hz, double bandwidth_ref_frequency_hz,
                       double spectral_exponent, bool kernel_truncation);

  std::vector<Constraint::Result> Apply(SolutionSpan& solutions, double time,
                                        std::ostream* stat_stream) override;

  /**
   * Set the weights to be used during smoothing.
   * @param weights An array of n_antenna x n_channel_blocks.
   */
  void SetWeights(const std::vector<double>& weights) override {
    weights_ = weights;
  }

  void SetSolutionWeights(
      const std::vector<std::vector<double>>& solution_weights) final {
    solution_weights_ = solution_weights;
  }

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  /**
   * Set antenna smoothness factors that can control the amount of smoothing
   * per antenna. One option is e.g. to let this depend on the distance to
   * the array centre. @note Because distance is the common use-case for these
   * factors, the factors work opposite of the dd smoothing factors: higher
   * values cause less smoothing.
   * @param antenna_factors vector where each element is a smoothing correction
   * factor for each antenna. A higher correction factor will perform stronger
   * smoothing in frequency direction.
   */
  void SetAntennaFactors(std::vector<double>&& antenna_factors) {
    antenna_factors_ = std::move(antenna_factors);
  }

  const std::vector<double>& GetAntennaFactors() const {
    return antenna_factors_;
  }

  /**
   * Sets an extra (cummulative) smoothing factor for directions (and dd
   * solution intervals). @note that these factors have opposite meaning
   * of those for @ref SetAntennaFactors().
   * @param dd_smoothing_factors should be of size NSolutions(), i.e.
   * it should have a factor for each direction and the possible
   * subsolutions per direction. The otherwise specified smoothing kernel size
   * is multiplied with these factors, which means that higher values will cause
   * stronger smoothing. If empty, all factors are assumed to be one.
   */
  void SetDdSmoothingFactors(std::vector<double> dd_smoothing_factors);

  const std::vector<double>& GetDdSmoothingFactors() const {
    return dd_smoothing_factors_;
  }

 private:
  using Smoother = KernelSmoother<std::complex<double>, double>;

  struct FitData {
    FitData(const std::vector<double>& frequencies,
            Smoother::KernelType kernel_type, double kernel_bandwidth,
            double bandwidth_ref_frequency_hz, double spectral_exponent,
            bool kernel_truncation)
        : smoother(frequencies, kernel_type, kernel_bandwidth,
                   bandwidth_ref_frequency_hz, spectral_exponent,
                   kernel_truncation),
          data(frequencies.size()),
          weight(frequencies.size()) {}

    Smoother smoother;
    std::vector<std::complex<double>> data;
    std::vector<double> weight;
  };
  std::vector<FitData> fit_data_;
  std::vector<double> frequencies_;
  std::vector<double> antenna_factors_;
  std::vector<double> dd_smoothing_factors_;
  std::vector<double> weights_;
  /// A weight array per solution. If given, these override @ref weights_.
  std::vector<std::vector<double>> solution_weights_;
  Smoother::KernelType kernel_type_;
  double bandwidth_;
  double bandwidth_ref_frequency_;
  double spectral_exponent_;
  bool kernel_truncation_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
