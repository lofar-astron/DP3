// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_DIAGONAL_ONLY_CONSTRAINT_H_
#define DP3_DDECAL_DIAGONAL_ONLY_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>

namespace dp3 {
namespace ddecal {

class DiagonalConstraint final : public Constraint {
 public:
  DiagonalConstraint(size_t pols_per_solution)
      : pols_per_solution_(pols_per_solution){};

  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions,
      [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) override {
    if (pols_per_solution_ == 4) {
      for (size_t ch = 0; ch < solutions.size(); ++ch) {
        for (size_t solIndex = 0; solIndex < solutions[ch].size();
             solIndex += 4) {
          solutions[ch][solIndex + 1] = 0.0;
          solutions[ch][solIndex + 2] = 0.0;
        }
      }
    }

    return std::vector<Constraint::Result>();
  }

 private:
  const size_t pols_per_solution_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
