// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "RotationAndDiagonalConstraint.h"
#include "RotationConstraint.h"

#include <cmath>

namespace {
bool dataIsValid(const std::complex<double>* const data,
                 const std::size_t size) {
  for (std::size_t i = 0; i < size; ++i) {
    if (std::isnan(data[i].real()) || std::isnan(data[i].imag())) {
      return false;
    }
  }
  return true;
}
}  // namespace

namespace dp3 {
namespace ddecal {

RotationAndDiagonalConstraint::RotationAndDiagonalConstraint()
    : _res(), _doRotationReference(false) {}

void RotationAndDiagonalConstraint::Initialize(
    size_t nAntennas, size_t nDirections,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, nDirections, frequencies);

  if (NDirections() != 1)  // TODO directions!
    throw std::runtime_error(
        "RotationAndDiagonalConstraint can't handle multiple directions yet");

  _res.resize(3);
  _res[0].vals.resize(NAntennas() * NChannelBlocks());
  _res[0].weights.resize(NAntennas() * NChannelBlocks());
  _res[0].axes = "ant,dir,freq";
  _res[0].dims.resize(3);
  _res[0].dims[0] = NAntennas();
  _res[0].dims[1] = NDirections();
  _res[0].dims[2] = NChannelBlocks();
  _res[0].name = "rotation";

  _res[1].vals.resize(NAntennas() * NChannelBlocks() * 2);
  _res[1].weights.resize(NAntennas() * NChannelBlocks() * 2);
  _res[1].axes = "ant,dir,freq,pol";
  _res[1].dims.resize(4);
  _res[1].dims[0] = NAntennas();
  _res[1].dims[1] = NDirections();
  _res[1].dims[2] = NChannelBlocks();
  _res[1].dims[3] = 2;
  _res[1].name = "amplitude";

  _res[2] = _res[1];
  _res[2].name = "phase";
}

void RotationAndDiagonalConstraint::SetWeights(
    const std::vector<double>& weights) {
  _res[0].weights = weights;  // TODO should be nDir times

  // Duplicate weights for two polarizations
  _res[1].weights.resize(NAntennas() * NDirections() * NChannelBlocks() * 2);
  size_t indexInWeights = 0;
  for (auto weight : weights) {
    _res[1].weights[indexInWeights++] = weight;  // TODO directions!
    _res[1].weights[indexInWeights++] = weight;
  }

  _res[2].weights = _res[1].weights;  // TODO directions!
}

void RotationAndDiagonalConstraint::SetDoRotationReference(
    const bool doRotationReference) {
  _doRotationReference = doRotationReference;
}

std::vector<Constraint::Result> RotationAndDiagonalConstraint::Apply(
    std::vector<std::vector<std::complex<double>>>& solutions, double,
    std::ostream* statStream) {
  if (statStream) *statStream << "[";  // begin channel
  double angle0 = std::nan("");
  for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
    if (statStream) *statStream << "[";  // begin antenna

    // First iterate over all antennas to find mean amplitudes, needed for
    // maxratio constraint below
    double amean = 0.0;
    double bmean = 0.0;
    for (size_t ant = 0; ant < NAntennas(); ++ant) {
      std::complex<double>* data = &(solutions[ch][4 * ant]);

      // Skip this antenna if has no valid data.
      if (!dataIsValid(data, 4)) {
        continue;
      }

      // Compute rotation
      double angle = RotationConstraint::get_rotation(data);
      // Restrict angle between -pi/2 and pi/2
      // Add 2pi to make sure that fmod doesn't see negative numbers
      angle = std::fmod(angle + 3.5 * M_PI, M_PI) - 0.5 * M_PI;

      // Right multiply solution with inverse rotation,
      // save only the diagonal
      // Use sin(-phi) == -sin(phi)
      std::complex<double> a = data[0] * cos(angle) - data[1] * sin(angle);
      std::complex<double> b = data[3] * cos(angle) + data[2] * sin(angle);

      double absa = std::abs(a);
      if (isfinite(absa)) {
        amean += absa;
      }
      double absb = std::abs(b);
      if (isfinite(absb)) {
        bmean += absb;
      }
    }
    amean /= NAntennas();
    bmean /= NAntennas();

    // Now iterate again to do the actual constraining
    bool diverged = false;
    for (size_t ant = 0; ant < NAntennas(); ++ant) {
      std::complex<double>* data = &(solutions[ch][4 * ant]);

      // Skip this antenna if has no valid data.
      if (!dataIsValid(data, 4)) {
        continue;
      }

      // Compute rotation
      double angle = RotationConstraint::get_rotation(data);

      // Restrict angle between -pi/2 and pi/2
      // Add 2pi to make sure that fmod doesn't see negative numbers
      angle = fmod(angle + 3.5 * M_PI, M_PI) - 0.5 * M_PI;

      // Right multiply solution with inverse rotation,
      // save only the diagonal
      // Use sin(-phi) == -sin(phi)
      std::complex<double> a = data[0] * cos(angle) - data[1] * sin(angle);
      std::complex<double> b = data[3] * cos(angle) + data[2] * sin(angle);

      // Constrain amplitudes to 1/maxratio < amp < maxratio
      double maxratio = 5.0;
      if (amean > 0.0) {
        if (std::abs(a) / amean < 1.0 / maxratio ||
            std::abs(a) / amean > maxratio) {
          diverged = true;
        }
        do {
          a *= 1.2;
        } while (std::abs(a) / amean < 1.0 / maxratio);
        do {
          a /= 1.2;
        } while (std::abs(a) / amean > maxratio);
      }
      if (bmean > 0.0) {
        if (std::abs(b) / bmean < 1.0 / maxratio ||
            std::abs(b) / bmean > maxratio) {
          diverged = true;
        }
        do {
          b *= 1.2;
        } while (std::abs(b) / bmean < 1.0 / maxratio);
        do {
          b /= 1.2;
        } while (std::abs(b) / bmean > maxratio);
      }

      if (_doRotationReference) {
        // Use the first station with a non-NaN angle as reference station
        // (for every chanblock), to work around unitary ambiguity
        if (std::isnan(angle0)) {
          angle0 = angle;
          angle = 0.;
        } else {
          angle -= angle0;
          angle = std::fmod(angle + 3.5 * M_PI, M_PI) - 0.5 * M_PI;
        }
      }

      _res[0].vals[ant * NChannelBlocks() + ch] = angle;  // TODO directions!
      _res[1].vals[ant * NChannelBlocks() * 2 + 2 * ch] = std::abs(a);
      _res[1].vals[ant * NChannelBlocks() * 2 + 2 * ch + 1] = std::abs(b);
      _res[2].vals[ant * NChannelBlocks() * 2 + 2 * ch] = std::arg(a);
      _res[2].vals[ant * NChannelBlocks() * 2 + 2 * ch + 1] = std::arg(b);

      // Do the actual constraining
      data[0] = a * cos(angle);
      data[1] = -a * sin(angle);
      data[2] = b * sin(angle);
      data[3] = b * cos(angle);
      if (statStream)
        *statStream << "[" << a.real() << "+" << a.imag() << "j," << b.real()
                    << "+" << b.imag() << "j," << angle << "]";
      // if (pd) cout<<angle;
      if (statStream && ant < NAntennas() - 1) *statStream << ",";
    }
    if (statStream) *statStream << "]";  // end antenna
    if (statStream && ch < NChannelBlocks() - 1) *statStream << ",";

    // If the maxratio constraint above was enforced for any antenna, set
    // weights of all antennas to a negative value for flagging later if desired
    if (diverged) {
      for (unsigned int ant = 0; ant < NAntennas(); ++ant) {
        _res[0].weights[ant * NChannelBlocks() + ch] = -1.0;
        _res[1].weights[ant * NChannelBlocks() * 2 + 2 * ch] = -1.0;
        _res[1].weights[ant * NChannelBlocks() * 2 + 2 * ch + 1] = -1.0;
        _res[2].weights[ant * NChannelBlocks() * 2 + 2 * ch] = -1.0;
        _res[2].weights[ant * NChannelBlocks() * 2 + 2 * ch + 1] = -1.0;
      }
    }
  }
  if (statStream) *statStream << "]\t";  // end channel

  return _res;
}

}  // namespace ddecal
}  // namespace dp3
