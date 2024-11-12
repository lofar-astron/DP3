// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_
#define DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_

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
                       double spectral_exponent);

  std::vector<Constraint::Result> Apply(SolutionSpan& solutions, double time,
                                        std::ostream* stat_stream) override;

  /**
   * Set the weights to be used during smoothing.
   * @param weights An array of n_antenna x n_channel_blocks.
   */
  void SetWeights(const std::vector<double>& weights) override {
    weights_ = weights;
  }

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  /**
   * Should be called after calling @ref Initialize().
   * @param antennaDistances vector where each element is a distance correction
   * factor for each antenna. A higher correction factor will perform stronger
   * smoothing in frequency direction.
   */
  void SetDistanceFactors(std::vector<double>&& antennaDistanceFactors);

  const std::vector<double>& GetDistanceFactors() const {
    return antenna_distance_factors_;
  }

 private:
  using Smoother = KernelSmoother<std::complex<double>, double>;

  struct FitData {
    FitData(const std::vector<double>& frequencies,
            Smoother::KernelType kernel_type, double kernel_bandwidth,
            double bandwidth_ref_frequency_hz, double spectral_exponent)
        : smoother(frequencies, kernel_type, kernel_bandwidth,
                   bandwidth_ref_frequency_hz, spectral_exponent),
          data(frequencies.size()),
          weight(frequencies.size()) {}

    Smoother smoother;
    std::vector<std::complex<double>> data;
    std::vector<double> weight;
  };
  std::vector<FitData> fit_data_;
  std::vector<double> frequencies_;
  std::vector<double> antenna_distance_factors_;
  std::vector<double> weights_;
  Smoother::KernelType kernel_type_;
  double bandwidth_;
  double bandwidth_ref_frequency_;
  double spectral_exponent_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
