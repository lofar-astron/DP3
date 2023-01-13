// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_PHASE_ONLY_CONSTRAINT_H_
#define DP3_DDECAL_PHASE_ONLY_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>

namespace dp3 {
namespace ddecal {

/**
 * @brief This class constrains the amplitudes of the solution to be unity, but
 * keeps the phase.
 */
class PhaseOnlyConstraint final : public Constraint {
 public:
  PhaseOnlyConstraint(){};

  std::vector<Constraint::Result> Apply(
      SolutionsSpan& solutions, [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) override {
    solutions /= xt::abs(solutions);

    return std::vector<Constraint::Result>();
  }
};

}  // namespace ddecal
}  // namespace dp3

#endif
