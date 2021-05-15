// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ROTATION_CONSTRAINT_H
#define ROTATION_CONSTRAINT_H

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

class RotationConstraint : public Constraint {
 public:
  RotationConstraint(){};

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex> >& solutions, double time,
      std::ostream* statStream);

  void Initialize(size_t nAntennas, size_t nDirections,
                  const std::vector<double>& frequencies) override;

  virtual void SetWeights(const std::vector<double>& weights);

  /* Compute the rotation from a 2x2 full jones solution */
  static double get_rotation(std::complex<double>* data);

 private:
  std::vector<Constraint::Result> _res;
};

}  // namespace ddecal
}  // namespace dp3

#endif
