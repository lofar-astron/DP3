// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Constraint.h"
#include "KernelSmoother.h"

#include <aocommon/parallelfor.h>

#ifndef SMOOTHNESS_CONSTRAINT_H
#define SMOOTHNESS_CONSTRAINT_H

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
  SmoothnessConstraint(double bandwidthHz, double bandwidthRefFrequencyHz);

  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double,
      std::ostream* statStream) final override;

  void SetWeights(const std::vector<double>& weights) final override {
    _weights = weights;
  }

  void Initialize(size_t nAntennas, size_t nDirections,
                  const std::vector<double>& frequencies) final override;

  /**
   * Should be called after calling @ref Initialize().
   * @param antennaDistances vector where each element is a distance correction
   * factor for each antenna. A higher correction factor will perform stronger
   * smoothing in frequency direction.
   */
  void SetDistanceFactors(std::vector<double>&& antennaDistanceFactors);

  const std::vector<double>& GetDistanceFactors() const {
    return _antennaDistanceFactors;
  }

  struct FitData {
    FitData(const std::vector<double>& frequencies,
            Smoother::KernelType kernelType, double kernelBandwidth,
            double bandwidthRefFrequencyHz)
        : smoother(frequencies, kernelType, kernelBandwidth,
                   bandwidthRefFrequencyHz),
          data(frequencies.size()),
          weight(frequencies.size()) {}

    Smoother smoother;
    std::vector<dcomplex> data;
    std::vector<double> weight;
  };
  std::vector<FitData> _fitData;
  std::vector<double> _frequencies;
  std::vector<double> _antennaDistanceFactors;
  std::vector<double> _weights;
  Smoother::KernelType _kernelType;
  double _bandwidth;
  double _bandwidthRefFrequencyHz;
  std::unique_ptr<aocommon::ParallelFor<size_t>> _loop;
};

}  // namespace ddecal
}  // namespace dp3

#endif
