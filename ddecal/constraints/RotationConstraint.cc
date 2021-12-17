// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "RotationConstraint.h"

#include <cmath>
#include <assert.h>

namespace dp3 {
namespace ddecal {

void RotationConstraint::Initialize(
    size_t nAntennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, solutions_per_direction, frequencies);

  if (NDirections() != 1)  // TODO directions!
    throw std::runtime_error(
        "RotationConstraint can't handle multiple directions yet");

  _res.resize(1);
  _res[0].vals.resize(NAntennas() * NChannelBlocks());
  _res[0].axes = "ant,dir,freq";
  _res[0].dims.resize(3);
  _res[0].dims[0] = NAntennas();
  _res[0].dims[1] = NDirections();
  _res[0].dims[2] = NChannelBlocks();
  _res[0].name = "rotation";
}

void RotationConstraint::SetWeights(const std::vector<double>& weights) {
  _res[0].weights = weights;  // TODO directions!
}

double RotationConstraint::GetRotation(std::complex<double>* data) {
  // Convert to circular
  dcomplex i(0, 1.);

  dcomplex ll = data[0] + data[3] - i * data[1] + i * data[2];
  dcomplex rr = data[0] + data[3] + i * data[1] - i * data[2];
  double angle = 0.5 * (arg(ll) - arg(rr));

  return angle;
}

std::vector<Constraint::Result> RotationConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double,
    std::ostream* /*statStream*/) {
  for (unsigned int ch = 0; ch < NChannelBlocks(); ++ch) {
    for (unsigned int ant = 0; ant < NAntennas(); ++ant) {
      // Compute rotation
      dcomplex* data = &(solutions[ch][4 * ant]);
      double angle = GetRotation(data);
      _res[0].vals[ant * NChannelBlocks() + ch] = angle;  // TODO directions!

      // Constrain the data
      data[0] = cos(angle);
      data[1] = -sin(angle);
      data[2] = -data[1];
      data[3] = data[0];
    }
  }

  return _res;
}

}  // namespace ddecal
}  // namespace dp3
