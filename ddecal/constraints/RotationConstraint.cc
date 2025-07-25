// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "RotationConstraint.h"

#include <cassert>
#include <cmath>

namespace dp3 {
namespace ddecal {

void RotationConstraint::Initialize(
    size_t nAntennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, solutions_per_direction, frequencies);

  Result& result = results_.emplace_back();
  result.vals.resize(NAntennas() * NSubSolutions() * NChannelBlocks());
  result.axes = "ant,dir,freq";
  result.dims.resize(3);
  result.dims[0] = NAntennas();
  result.dims[1] = NSubSolutions();
  result.dims[2] = NChannelBlocks();
  result.name = "rotation";
  // This constraint supports dd solution intervals, but the hdf5 writer
  // code does not yet support it for constraint results.
  if (NSubSolutions() != NDirections()) {
    throw std::runtime_error(
        "The rotation constraint does not support direction-dependent "
        "intervals");
  }
}

void RotationConstraint::SetWeights(const std::vector<double>& weights) {
  Result& result = results_.front();
  result.weights.clear();
  for (size_t i = 0; i != NSubSolutions(); ++i) {
    result.weights.insert(result.weights.end(), weights.begin(), weights.end());
  }
}

std::vector<Constraint::Result> RotationConstraint::Apply(
    SolutionSpan& solutions, double,
    [[maybe_unused]] std::ostream* statStream) {
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) == 4);  // 2x2 full jones solutions
  Result& result = results_.front();
  for (size_t sub_solution = 0; sub_solution != NSubSolutions();
       ++sub_solution) {
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      for (size_t ant = 0; ant < NAntennas(); ++ant) {
        // Compute rotation
        std::complex<double>* data = &solutions(ch, ant, sub_solution, 0);
        const double angle = FitRotation(data);
        const size_t index =
            ((ant * NSubSolutions()) + sub_solution) * NChannelBlocks() + ch;
        result.vals[index] = angle;

        // Constrain the data
        RotationConstraint::SetRotation(data, angle);
      }
    }
  }

  return results_;
}

}  // namespace ddecal
}  // namespace dp3
