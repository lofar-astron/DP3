// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "KernelSmoother.h"
#include "SmoothnessConstraint.h"
#include <aocommon/parallelfor.h>
#include <boost/make_unique.hpp>

SmoothnessConstraint::SmoothnessConstraint(double bandwidthHz,
                                           double bandwidthRefFrequencyHz)
    : _kernelType(Smoother::GaussianKernel),
      _bandwidth(bandwidthHz),
      _bandwidthRefFrequencyHz(bandwidthRefFrequencyHz) {}

void SmoothnessConstraint::Initialize(size_t nAntennas, size_t nDirections,
                                      const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, nDirections, frequencies);
  _frequencies = frequencies;
}

void SmoothnessConstraint::SetDistanceFactors(
    std::vector<double>&& antennaDistanceFactors) {
  _antennaDistanceFactors = std::move(antennaDistanceFactors);
  if (!_loop) {
    _loop = boost::make_unique<aocommon::ParallelFor<size_t>>(NThreads());
  }
  for (size_t i = 0; i != NThreads(); ++i)
    _fitData.emplace_back(_frequencies, _kernelType, _bandwidth,
                          _bandwidthRefFrequencyHz);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double, std::ostream*) {
  const size_t nPol = solutions.front().size() / (NAntennas() * NDirections());

  _loop->Run(
      0, NAntennas() * NDirections(), [&](size_t antDirIndex, size_t thread) {
        size_t antIndex = antDirIndex / NDirections();
        for (size_t pol = 0; pol != nPol; ++pol) {
          size_t solutionIndex = antDirIndex * nPol + pol;
          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            // Flag channels where calibration yielded inf or nan
            if (isfinite(solutions[ch][solutionIndex])) {
              _fitData[thread].data[ch] = solutions[ch][solutionIndex];
              _fitData[thread].weight[ch] =
                  _weights[antIndex * NChannelBlocks() + ch];
            } else {
              _fitData[thread].data[ch] = 0.0;
              _fitData[thread].weight[ch] = 0.0;
            }
          }

          _fitData[thread].smoother.Smooth(_fitData[thread].data.data(),
                                           _fitData[thread].weight.data(),
                                           _antennaDistanceFactors[antIndex]);

          for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
            solutions[ch][solutionIndex] = _fitData[thread].data[ch];
          }
        }
      });

  return std::vector<Constraint::Result>();
}
