// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SCALAR_SOLVER_H
#define DDECAL_SCALAR_SOLVER_H

#include "Constraint.h"
#include "SolverBase.h"
#include "SolverBuffer.h"
#include "Stopwatch.h"

#include <complex>
#include <vector>
#include <memory>

namespace DP3 {
namespace DPPP {

class ScalarSolver : public SolverBase {
 public:
  ScalarSolver() : SolverBase() {}

  SolveResult Solve(const std::vector<Complex*>&, const std::vector<float*>&,
                    const std::vector<std::vector<Complex*>>&,
                    std::vector<std::vector<DComplex>>&, double,
                    std::ostream*) override;

 private:
  void PerformIteration(size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
                        std::vector<Matrix>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions);
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DDECAL_SCALAR_SOLVER_H
