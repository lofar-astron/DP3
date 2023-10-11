// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Constraint.h"
#include "KernelSmoother.h"

#include <aocommon/dynamicfor.h>

#ifndef DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_
#define DP3_DDECAL_SMOOTHNESS_CONSTRAINT_H_

namespace dp3 {
namespace ddecal {

class SmoothnessConstraint final : public Constraint {
 public:
  typedef std::complex<double> dcomplex;
  typedef KernelSmoother<dcomplex, double> Smoother;

  /**
   * @param bandwidth_hz Size of the kernel (smoothing strength)
   * @param bandwidth_ref_frequency_hz Reference frequency for the kernel size,
   * may be zero to have a constant kernel size over frequency.
   * More details on this are in the KernelSmoother documentation.
   */
  SmoothnessConstraint(double bandwidth_hz, double bandwidth_ref_frequency_hz);

  std::vector<Constraint::Result> Apply(SolutionSpan& solutions, double time,
                                        std::ostream* stat_stream) override;

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
};

}  // namespace ddecal
}  // namespace dp3

#endif
