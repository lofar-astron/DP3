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

#include "Constraint.h"
#include "KernelSmoother.h"

#include <aocommon/parallelfor.h>

#ifndef SMOOTHNESS_CONSTRAINT_H
#define SMOOTHNESS_CONSTRAINT_H

class SmoothnessConstraint : public Constraint {
 public:
  typedef std::complex<double> dcomplex;
  typedef KernelSmoother<dcomplex, double> Smoother;

  SmoothnessConstraint(double bandwidthHz);

  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double,
      std::ostream* statStream) final override;

  void SetWeights(const std::vector<double>& weights) final override {
    _weights = weights;
  }

  void Initialize(const double* frequencies);

  virtual void InitializeDimensions(size_t nAntennas, size_t nDirections,
                                    size_t nChannelBlocks) final override;

  struct FitData {
    FitData(const double* frequencies, size_t n,
            Smoother::KernelType kernelType, double kernelBandwidth)
        : smoother(frequencies, n, kernelType, kernelBandwidth),
          data(n),
          weight(n) {}

    Smoother smoother;
    std::vector<dcomplex> data;
    std::vector<double> weight;
  };
  std::vector<FitData> _fitData;
  std::vector<double> _frequencies, _weights;
  Smoother::KernelType _kernelType;
  double _bandwidth;
  std::unique_ptr<aocommon::ParallelFor<size_t>> _loop;
};

#endif
