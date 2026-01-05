// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_
#define DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_

#include "Constraint.h"

#include "base/CalType.h"

#include <vector>
#include <optional>
#include <ostream>

namespace dp3 {
namespace ddecal {

class RotationAndDiagonalConstraint final : public Constraint {
 public:
  RotationAndDiagonalConstraint(base::CalType diagonal_solution_type);

  std::vector<ConstraintResult> Apply(SolutionSpan& solutions, double time,
                                      std::ostream* statStream) override;

  void Initialize(size_t n_antennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  void SetWeights(const std::vector<double>& weights) override;

  void SetDoRotationReference(bool do_rotation_reference);

 private:
  void FitRotationAndDiagonal(
      const std::complex<double>* data, double& angle,
      std::array<std::complex<double>, 2>& diagonal) const;

  template <size_t PolCount>
  void SetChannelWeights(std::vector<double>& values, size_t channel,
                         double new_value) const;

  std::vector<ConstraintResult> results_;
  bool do_rotation_reference_;
  base::CalType diagonal_solution_type_;
};

/**
 * Constrain the diagonal in accordance with the supplied mode. The mode
 * should be a diagonal, scalar or rotational mode. In case
 * of rotational, the diagonal is set to unity.
 */
void ConstrainDiagonal(std::array<std::complex<double>, 2>& diagonal,
                       base::CalType mode);

std::vector<ConstraintResult> MakeDiagonalResults(
    size_t n_antennas, size_t n_sub_solutions, size_t n_channels,
    base::CalType diagonal_solution_type);

inline void StoreDiagonal(ConstraintResult* results,
                          const std::array<std::complex<double>, 2>& diagonal,
                          size_t channel, size_t antenna, size_t subsolution,
                          size_t n_channels, size_t n_sub_solutions,
                          base::CalType diagonal_solution_type) {
  // This is a bit more verbose than necessary, but using a single
  // switch statement instead of multiple conditionals is probably
  // most efficient.
  using base::CalType;
  const size_t index_part =
      (antenna * n_sub_solutions + subsolution) * n_channels + channel;
  switch (diagonal_solution_type) {
    case CalType::kDiagonal:
      results[0].vals[index_part * 2] = std::abs(diagonal[0]);
      results[0].vals[index_part * 2 + 1] = std::abs(diagonal[1]);
      results[1].vals[index_part * 2] = std::arg(diagonal[0]);
      results[1].vals[index_part * 2 + 1] = std::arg(diagonal[1]);
      break;
    case CalType::kDiagonalAmplitude:
      results[0].vals[index_part * 2] = std::abs(diagonal[0]);
      results[0].vals[index_part * 2 + 1] = std::abs(diagonal[1]);
      break;
    case CalType::kDiagonalPhase:
      results[0].vals[index_part * 2] = std::arg(diagonal[0]);
      results[0].vals[index_part * 2 + 1] = std::arg(diagonal[1]);
      break;
    case CalType::kScalar:
      results[0].vals[index_part] = std::abs(diagonal[0]);
      results[1].vals[index_part] = std::arg(diagonal[0]);
      break;
    case CalType::kScalarAmplitude:
      results[0].vals[index_part] = std::abs(diagonal[0]);
      break;
    case CalType::kScalarPhase:
      results[0].vals[index_part] = std::arg(diagonal[0]);
      break;
    default:
      assert(false);
  }
}

}  // namespace ddecal
}  // namespace dp3

#endif
