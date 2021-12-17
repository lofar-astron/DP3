// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ROTATION_CONSTRAINT_H_
#define DP3_DDECAL_ROTATION_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

class RotationConstraint final : public Constraint {
 public:
  RotationConstraint(){};

  std::vector<Result> Apply(std::vector<std::vector<dcomplex>>& solutions,
                            double time, std::ostream* statStream) override;

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  void SetWeights(const std::vector<double>& weights) override;

  /* Compute the rotation from a 2x2 full jones solution */
  static double GetRotation(std::complex<double>* data);

 private:
  std::vector<Constraint::Result> _res;
};

}  // namespace ddecal
}  // namespace dp3

#endif
