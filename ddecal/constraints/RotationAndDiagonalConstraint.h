// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ROTATIONANDDIAGONAL_CONSTRAINT_H
#define ROTATIONANDDIAGONAL_CONSTRAINT_H

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace dp3 {

class RotationAndDiagonalConstraint : public Constraint {
 public:
  RotationAndDiagonalConstraint();

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex> >& solutions, double time,
      std::ostream* statStream);

  virtual void InitializeDimensions(size_t nAntennas, size_t nDirections,
                                    size_t nChannelBlocks);

  virtual void SetWeights(const std::vector<double>& weights);

  void SetDoRotationReference(bool doRotationReference);

 private:
  std::vector<Constraint::Result> _res;
  bool _doRotationReference;
};

}  // namespace dp3

#endif
