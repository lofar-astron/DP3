// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_AMPLITUDE_ONLY_CONSTRAINT_H_
#define DP3_DDECAL_AMPLITUDE_ONLY_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>

namespace dp3 {
namespace ddecal {

/**
 * @brief This class constrains the phases of the solution to be zero, but
 * keeps the amplitude information.
 */
class AmplitudeOnlyConstraint final : public Constraint {
 public:
  std::vector<Constraint::Result> Apply(
      SolutionSpan& solutions, double,
      [[maybe_unused]] std::ostream* stat_stream) override {
    solutions = xt::abs(solutions);

    return std::vector<Constraint::Result>();
  }
};

}  // namespace ddecal
}  // namespace dp3

#endif
