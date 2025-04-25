// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "RotationAndDiagonalConstraint.h"
#include "RotationConstraint.h"

#include <cmath>

using dp3::base::CalType;

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

/// Adjusts the magnitude of value to fall within a specified range from
/// mean_abs. Returns true if the initial value was already within the range,
/// false otherwise.
bool Limit(std::complex<double>& value, double mean_abs, double max_ratio) {
  if (mean_abs > 0.0) {
    if (std::abs(value) / mean_abs < 1.0 / max_ratio ||
        std::abs(value) / mean_abs > max_ratio) {
      return true;
    }
    do {
      value *= 1.2;
    } while (std::abs(value) / mean_abs < 1.0 / max_ratio);
    do {
      value /= 1.2;
    } while (std::abs(value) / mean_abs > max_ratio);
  }
  return false;
}

}  // namespace

namespace dp3 {
namespace ddecal {

RotationAndDiagonalConstraint::RotationAndDiagonalConstraint(
    CalType diagonal_solution_type)
    : results_(),
      do_rotation_reference_(false),
      diagonal_solution_type_(diagonal_solution_type) {
  if (diagonal_solution_type != CalType::kDiagonal &&
      diagonal_solution_type != CalType::kDiagonalAmplitude &&
      diagonal_solution_type != CalType::kDiagonalPhase &&
      diagonal_solution_type != CalType::kScalar &&
      diagonal_solution_type != CalType::kScalarAmplitude &&
      diagonal_solution_type != CalType::kScalarPhase) {
    throw std::invalid_argument(
        "The diagonal mode of the rotation-and-diagonal-constraint must be "
        "diagonal, diagonalamplitude, diagonalphase, scalar, scalaramplitude "
        "or scalarphase");
  }
}

void ConstrainDiagonal(std::array<std::complex<double>, 2>& diagonal,
                       CalType mode) {
  switch (mode) {
    case CalType::kDiagonal:
      break;
    case CalType::kDiagonalAmplitude:
      diagonal[0] = std::abs(diagonal[0]);
      diagonal[1] = std::abs(diagonal[1]);
      break;
    case CalType::kDiagonalPhase:
      diagonal[0] = diagonal[0] / std::abs(diagonal[0]);
      diagonal[1] = diagonal[1] / std::abs(diagonal[1]);
      break;
    case CalType::kScalar:
      diagonal[0] = (diagonal[0] + diagonal[1]) * 0.5;
      diagonal[1] = diagonal[0];
      break;
    case CalType::kScalarAmplitude:
      diagonal[0] = std::abs((diagonal[0] + diagonal[1]) * 0.5);
      diagonal[1] = diagonal[0];
      break;
    case CalType::kScalarPhase:
      diagonal[0] = (diagonal[0] + diagonal[1]) * 0.5;
      diagonal[0] = diagonal[0] / std::abs(diagonal[0]);
      diagonal[1] = diagonal[0];
      break;
    case CalType::kRotation:
      diagonal[0] = 1.0;
      diagonal[1] = 1.0;
      break;
    case CalType::kFullJones:
    case CalType::kRotationAndDiagonal:
    case CalType::kTec:
    case CalType::kTecAndPhase:
    case CalType::kTecScreen:
      assert(false);
      break;
  }
}

std::vector<Constraint::Result> MakeDiagonalResults(
    size_t n_antennas, size_t n_sub_solutions, size_t n_channels,
    CalType diagonal_solution_type) {
  const bool is_scalar = diagonal_solution_type == CalType::kScalar ||
                         diagonal_solution_type == CalType::kScalarAmplitude ||
                         diagonal_solution_type == CalType::kScalarPhase;
  const size_t n_diagonal_parameters = is_scalar ? 1 : 2;
  Constraint::Result template_result;
  template_result.vals.resize(n_antennas * n_sub_solutions * n_channels *
                              n_diagonal_parameters);
  template_result.weights.resize(n_antennas * n_sub_solutions * n_channels *
                                 n_diagonal_parameters);
  if (n_diagonal_parameters != 1) {
    template_result.axes = "ant,dir,freq,pol";
    template_result.dims.resize(4);
    template_result.dims[3] = n_diagonal_parameters;
  } else {
    template_result.axes = "ant,dir,freq";
    template_result.dims.resize(3);
  }
  template_result.dims[0] = n_antennas;
  template_result.dims[1] = n_sub_solutions;
  template_result.dims[2] = n_channels;

  std::vector<Constraint::Result> result;
  const bool has_amplitude =
      diagonal_solution_type == CalType::kDiagonal ||
      diagonal_solution_type == CalType::kDiagonalAmplitude ||
      diagonal_solution_type == CalType::kScalar ||
      diagonal_solution_type == CalType::kScalarAmplitude;
  if (has_amplitude) {
    Constraint::Result& amplitude_result = result.emplace_back();
    amplitude_result = template_result;
    amplitude_result.name = "amplitude";
  }

  const bool has_phase = diagonal_solution_type == CalType::kDiagonal ||
                         diagonal_solution_type == CalType::kDiagonalPhase ||
                         diagonal_solution_type == CalType::kScalar ||
                         diagonal_solution_type == CalType::kScalarPhase;
  if (has_phase) {
    Constraint::Result& phase_result = result.emplace_back();
    // The template is no longer necessary, so use its allocation (by moving)
    phase_result = std::move(template_result);
    phase_result.name = "phase";
  }
  return result;
}

void RotationAndDiagonalConstraint::Initialize(
    size_t nAntennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, solutions_per_direction, frequencies);

  // This constraint supports dd solution intervals, but the hdf5 writer
  // code does not yet support it for constraint results.
  if (NSubSolutions() != NDirections()) {
    throw std::runtime_error(
        "The rotation-and-diagonal constraint does not support "
        "direction-dependent "
        "intervals");
  }

  Result& rotation_result = results_.emplace_back();
  rotation_result.vals.resize(NAntennas() * NSubSolutions() * NChannelBlocks());
  rotation_result.weights.resize(NAntennas() * NSubSolutions() *
                                 NChannelBlocks());
  rotation_result.axes = "ant,dir,freq";
  rotation_result.dims.resize(3);
  rotation_result.dims[0] = NAntennas();
  rotation_result.dims[1] = NSubSolutions();
  rotation_result.dims[2] = NChannelBlocks();
  rotation_result.name = "rotation";

  std::vector<Result> diagonal_result = MakeDiagonalResults(
      NAntennas(), NSubSolutions(), NChannelBlocks(), diagonal_solution_type_);
  for (Result& result : diagonal_result)
    results_.emplace_back(std::move(result));
}

void RotationAndDiagonalConstraint::SetWeights(
    const std::vector<double>& weights) {
  // weights is nAntennas * nChannelBlocks
  results_[0].weights = weights;  // TODO should be nInterval times

  // Duplicate weights for one or two polarizations
  const size_t n_polarizations = GetNPolarizations(diagonal_solution_type_);
  results_[1].weights.resize(weights.size() * n_polarizations);
  size_t index_in_weights = 0;
  for (double weight : weights) {
    for (size_t p = 0; p != n_polarizations; ++p) {
      results_[1].weights[index_in_weights + p] =
          weight;  // TODO directions / intervals!
    }
    index_in_weights += n_polarizations;
  }

  if (results_.size() > 2)
    results_[2].weights = results_[1].weights;  // TODO directions / intervals!
}

void RotationAndDiagonalConstraint::SetDoRotationReference(
    const bool doRotationReference) {
  do_rotation_reference_ = doRotationReference;
}

void RotationAndDiagonalConstraint::FitRotationAndDiagonal(
    const std::complex<double>* data, double& angle,
    std::array<std::complex<double>, 2>& diagonal) const {
  // Compute rotation
  angle = RotationConstraint::FitRotation(data);
  // Restrict angle between -pi/2 and pi/2
  // Add 2pi to make sure that fmod doesn't see negative numbers
  angle = std::fmod(angle + 3.5 * M_PI, M_PI) - 0.5 * M_PI;

  // Right multiply solution with inverse rotation,
  // save only the diagonal
  // Use sin(-phi) == -sin(phi)
  diagonal[0] = data[0] * std::cos(angle) - data[1] * std::sin(angle);
  diagonal[1] = data[3] * std::cos(angle) + data[2] * std::sin(angle);
}

template <size_t PolCount>
void RotationAndDiagonalConstraint::SetChannelWeights(
    std::vector<double>& values, size_t channel, double new_value) const {
  for (size_t antenna = 0; antenna != NAntennas(); ++antenna) {
    for (size_t p = 0; p != PolCount; ++p) {
      values[(antenna * NChannelBlocks() + channel) * PolCount + p] = new_value;
    }
  }
}

std::vector<Constraint::Result> RotationAndDiagonalConstraint::Apply(
    SolutionSpan& solutions, double,
    [[maybe_unused]] std::ostream* statStream) {
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) == 4);  // 2x2 full jones solutions.

  for (size_t sub_solution = 0; sub_solution != NSubSolutions();
       ++sub_solution) {
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      // First iterate over all antennas to find mean amplitudes, needed for
      // maxratio constraint below
      double a_mean = 0.0;
      double b_mean = 0.0;
      for (size_t ant = 0; ant < NAntennas(); ++ant) {
        std::complex<double>* data = &(solutions(ch, ant, sub_solution, 0));

        // Skip this antenna if has no valid data.
        if (!dataIsValid(data, 4)) {
          continue;
        }

        double angle;
        std::array<std::complex<double>, 2> diagonal;
        FitRotationAndDiagonal(data, angle, diagonal);

        const double abs_a = std::abs(diagonal[0]);
        if (std::isfinite(abs_a)) {
          a_mean += abs_a;
        }
        const double abs_b = std::abs(diagonal[1]);
        if (std::isfinite(abs_b)) {
          b_mean += abs_b;
        }
      }
      a_mean /= NAntennas();
      b_mean /= NAntennas();

      double angle0 = std::numeric_limits<double>::quiet_NaN();

      // Now iterate again to do the actual constraining
      bool diverged = false;
      for (size_t ant = 0; ant != NAntennas(); ++ant) {
        std::complex<double>* data = &solutions(ch, ant, sub_solution, 0);

        // Skip this antenna if has no valid data.
        if (!dataIsValid(data, 4)) {
          continue;
        }

        double angle;
        std::array<std::complex<double>, 2> diagonal;
        FitRotationAndDiagonal(data, angle, diagonal);

        // Constrain amplitudes to 1/maxratio < amp < maxratio
        double maxratio = 5.0;
        if (a_mean > 0.0 && Limit(diagonal[0], a_mean, maxratio)) {
          diverged = true;
        }
        if (b_mean > 0.0 && Limit(diagonal[1], a_mean, maxratio)) {
          diverged = true;
        }

        ConstrainDiagonal(diagonal, diagonal_solution_type_);

        if (do_rotation_reference_) {
          // Use the first station with a non-NaN angle as reference station
          // (for every chanblock), to work around unitary ambiguity
          if (std::isnan(angle0)) {
            angle0 = angle;
            angle = 0.0;
          } else {
            angle -= angle0;
            angle = std::fmod(angle + 3.5 * M_PI, M_PI) - 0.5 * M_PI;
          }
        }

        const size_t index =
            ((ant * NSubSolutions()) + sub_solution * NChannelBlocks()) + ch;
        results_[0].vals[index] = angle;
        StoreDiagonal(&results_[1], diagonal, ch, ant, sub_solution,
                      NChannelBlocks(), NSubSolutions(),
                      diagonal_solution_type_);

        // Do the actual constraining
        data[0] = diagonal[0] * std::cos(angle);
        data[1] = -diagonal[0] * std::sin(angle);
        data[2] = diagonal[1] * std::sin(angle);
        data[3] = diagonal[1] * std::cos(angle);
      }

      // If the maxratio constraint above was enforced for any antenna, set
      // weights of all antennas to a negative value for flagging later if
      // desired
      if (diverged) {
        SetChannelWeights<1>(results_[0].weights, ch, -1.0);
        for (size_t i = 1; i != results_.size(); ++i) {
          if (GetNPolarizations(diagonal_solution_type_) == 1)
            SetChannelWeights<1>(results_[i].weights, ch, -1.0);
          else
            SetChannelWeights<2>(results_[i].weights, ch, -1.0);
        }
      }
    }
  }

  return results_;
}

}  // namespace ddecal
}  // namespace dp3
