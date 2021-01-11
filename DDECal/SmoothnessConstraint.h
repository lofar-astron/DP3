// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Constraint.h"
#include "KernelSmoother.h"

#include <aocommon/parallelfor.h>

#ifndef SMOOTHNESS_CONSTRAINT_H
#define SMOOTHNESS_CONSTRAINT_H

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

  /**
   * Should be called after constructing.
   * @param frequencies list of channel frequencies in Hz
   * @param antennaDistances vector where each element is a distance correction
   * factor for each antenna. A higher correction factor will perform stronger
   * smoothing in frequency direction.
   */
  void Initialize(const double* frequencies,
                  std::vector<double> antennaDistanceFactors);

  virtual void InitializeDimensions(size_t nAntennas, size_t nDirections,
                                    size_t nChannelBlocks) final override;

  struct FitData {
    FitData(const double* frequencies, size_t n,
            Smoother::KernelType kernelType, double kernelBandwidth,
            double bandwidthRefFrequencyHz)
        : smoother(frequencies, n, kernelType, kernelBandwidth,
                   bandwidthRefFrequencyHz),
          data(n),
          weight(n) {}

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

#endif
