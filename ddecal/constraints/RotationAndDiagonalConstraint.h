// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_
#define DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_

#include "Constraint.h"

#include "../../base/CalType.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

void ConstrainDiagonal(std::array<std::complex<double>, 2>& diagonal,
                       base::CalType mode);

class RotationAndDiagonalConstraint final : public Constraint {
 public:
  RotationAndDiagonalConstraint(base::CalType diagonal_solution_type);

  std::vector<Result> Apply(SolutionSpan& solutions, double time,
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
  void StoreDiagonal(const std::array<std::complex<double>, 2>& diagonal,
                     size_t antenna, size_t channel);

  template <size_t PolCount>
  void SetChannelWeights(std::vector<double>& values, size_t channel,
                         double new_value) const;

  std::vector<Constraint::Result> results_;
  bool do_rotation_reference_;
  base::CalType diagonal_solution_type_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
