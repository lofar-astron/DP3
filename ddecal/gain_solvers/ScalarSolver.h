// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SCALAR_SOLVER_H
#define DDECAL_SCALAR_SOLVER_H

#include "RegularSolverBase.h"

namespace dp3 {
namespace ddecal {

class ScalarSolver final : public RegularSolverBase {
 public:
  ScalarSolver() : RegularSolverBase() {}

  SolveResult Solve(const SolverBuffer& solver_buffer,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 1; }

 private:
  void PerformIteration(const SolverBuffer& solver_buffer,
                        size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
                        std::vector<Matrix>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions,
                        double iterationfraction, double solverprecision);
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_SCALAR_SOLVER_H
