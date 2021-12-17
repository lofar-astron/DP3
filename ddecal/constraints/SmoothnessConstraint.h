// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Constraint.h"
#include "KernelSmoother.h"

#include <aocommon/parallelfor.h>

#ifndef DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_
#define DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_

namespace dp3 {
namespace ddecal {

class SmoothnessConstraint : public Constraint {
 public:
  typedef std::complex<double> dcomplex;
  typedef KernelSmoother<dcomplex, double> Smoother;

  /**
   * @param bandwidthRefFrequencyHz may be zero to have a constant kernel size
   * over frequency.
   */
  SmoothnessConstraint(double bandwidth_hz, double bandwidth_ref_frequency_hz);

  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* stat_stream) final override;

  void SetWeights(const std::vector<double>& weights) final override {
    weights_ = weights;
  }

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) final override;

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

  struct FitData {
    FitData(const std::vector<double>& frequencies,
            Smoother::KernelType kernel_type, double kernel_bandwidth,
            double bandwidth_ref_frequency_hz)
        : smoother(frequencies, kernel_type, kernel_bandwidth,
                   bandwidth_ref_frequency_hz),
          data(frequencies.size()),
          weight(frequencies.size()) {}

    Smoother smoother;
    std::vector<dcomplex> data;
    std::vector<double> weight;
  };
  std::vector<FitData> fit_data_;
  std::vector<double> frequencies_;
  std::vector<double> antenna_distance_factors_;
  std::vector<double> weights_;
  Smoother::KernelType kernel_type_;
  double bandwidth_;
  double bandwidth_ref_frequency_;
  std::unique_ptr<aocommon::ParallelFor<size_t>> loop_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
