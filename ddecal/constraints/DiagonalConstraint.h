// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_DIAGONAL_ONLY_CONSTRAINT_H_
#define DP3_DDECAL_DIAGONAL_ONLY_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>

#include <xtensor/xview.hpp>

namespace dp3 {
namespace ddecal {

class DiagonalConstraint final : public Constraint {
 public:
  DiagonalConstraint(size_t pols_per_solution)
      : pols_per_solution_(pols_per_solution){};

  std::vector<Constraint::Result> Apply(
      SolutionSpan& solutions, [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) override {
    if (pols_per_solution_ == 4) {
      xt::view(solutions, xt::all(), xt::all(), xt::all(), xt::range(1, 3)) =
          0.0;
    }

    return std::vector<Constraint::Result>();
  }

 private:
  const size_t pols_per_solution_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
