// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_DIAGONAL_SOLVER_H
#define DDECAL_DIAGONAL_SOLVER_H

#include "Constraint.h"
#include "SolverBase.h"

#include <complex>
#include <vector>
#include <memory>

namespace DP3 {
namespace DPPP {

class DiagonalSolver : public SolverBase {
 public:
  DiagonalSolver() : SolverBase() {}

  SolveResult Solve(
      const std::vector<Complex*>& unweighted_data,
      const std::vector<float*>& weights,
      const std::vector<std::vector<Complex*>>& unweighted_model_data,
      std::vector<std::vector<DComplex>>& solutions, double time,
      std::ostream* stat_stream) override;

 private:
  void PerformIteration(size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
                        std::vector<std::vector<Complex>>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions);
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DDECAL_DIAGONAL_SOLVER_H
