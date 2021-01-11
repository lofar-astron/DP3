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

void SmoothnessConstraint::Initialize(
    const double* frequencies, std::vector<double> antennaDistanceFactors) {
  _frequencies.assign(frequencies, frequencies + _nChannelBlocks);
  _antennaDistanceFactors = std::move(antennaDistanceFactors);
  if (!_loop) {
    _loop = boost::make_unique<aocommon::ParallelFor<size_t>>(_nThreads);
  }
  for (size_t i = 0; i != _nThreads; ++i)
    _fitData.emplace_back(_frequencies.data(), _frequencies.size(), _kernelType,
                          _bandwidth, _bandwidthRefFrequencyHz);
}

void SmoothnessConstraint::InitializeDimensions(size_t nAntennas,
                                                size_t nDirections,
                                                size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double, std::ostream*) {
  const size_t nPol = solutions.front().size() / (_nAntennas * _nDirections);

  _loop->Run(
      0, _nAntennas * _nDirections, [&](size_t antDirIndex, size_t thread) {
        size_t antIndex = antDirIndex / _nDirections;
        for (size_t pol = 0; pol != nPol; ++pol) {
          size_t solutionIndex = antDirIndex * nPol + pol;
          for (size_t ch = 0; ch != _nChannelBlocks; ++ch) {
            // Flag channels where calibration yielded inf or nan
            if (isfinite(solutions[ch][solutionIndex])) {
              _fitData[thread].data[ch] = solutions[ch][solutionIndex];
              _fitData[thread].weight[ch] =
                  _weights[antIndex * _nChannelBlocks + ch];
            } else {
              _fitData[thread].data[ch] = 0.0;
              _fitData[thread].weight[ch] = 0.0;
            }
          }

          _fitData[thread].smoother.Smooth(_fitData[thread].data.data(),
                                           _fitData[thread].weight.data(),
                                           _antennaDistanceFactors[antIndex]);

          for (size_t ch = 0; ch != _nChannelBlocks; ++ch) {
            solutions[ch][solutionIndex] = _fitData[thread].data[ch];
          }
        }
      });

  return std::vector<Constraint::Result>();
}
