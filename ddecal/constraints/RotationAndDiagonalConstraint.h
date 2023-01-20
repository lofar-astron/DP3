// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_
#define DP3_DDECAL_ROTATIONANDDIAGONAL_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

class RotationAndDiagonalConstraint final : public Constraint {
 public:
  RotationAndDiagonalConstraint();

  std::vector<Result> Apply(std::vector<std::vector<dcomplex>>& solutions,
                            double time, std::ostream* statStream) override;

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  void SetWeights(const std::vector<double>& weights) override;

  void SetDoRotationReference(bool doRotationReference);

 private:
  std::vector<Constraint::Result> _res;
  bool _doRotationReference;
};

}  // namespace ddecal
}  // namespace dp3

#endif
