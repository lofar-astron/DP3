// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "KernelSmoother.h"
#include "SmoothnessConstraint.h"
#include <aocommon/parallelfor.h>

SmoothnessConstraint::SmoothnessConstraint(double bandwidthHz)
    : _kernelType(Smoother::GaussianKernel), _bandwidth(bandwidthHz) {}

void SmoothnessConstraint::Initialize(const double* frequencies) {
  _frequencies.assign(frequencies, frequencies + _nChannelBlocks);
  if (!_loop) _loop.reset(new aocommon::ParallelFor<size_t>(_nThreads));
  for (size_t i = 0; i != _nThreads; ++i)
    _fitData.emplace_back(_frequencies.data(), _frequencies.size(), _kernelType,
                          _bandwidth);
}

void SmoothnessConstraint::InitializeDimensions(size_t nAntennas,
                                                size_t nDirections,
                                                size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double, std::ostream*) {
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
                                           _fitData[thread].weight.data());

          for (size_t ch = 0; ch != _nChannelBlocks; ++ch) {
            solutions[ch][solutionIndex] = _fitData[thread].data[ch];
          }
        }
      });

  return std::vector<Constraint::Result>();
}
