// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SCALAR_SOLVER_H
#define DDECAL_SCALAR_SOLVER_H

#include "SolverBase.h"

#include <complex>
#include <vector>

namespace DP3 {
namespace DPPP {

class ScalarSolver final : public SolverBase {
 public:
  ScalarSolver() : SolverBase() {}

  SolveResult Solve(const std::vector<Complex*>& unweighted_data,
                    const std::vector<float*>& weights,
                    std::vector<std::vector<Complex*>>&& unweighted_model_data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

 private:
  void PerformIteration(size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
                        std::vector<Matrix>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions,
                        double iterationfraction, double solverprecision);
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DDECAL_SCALAR_SOLVER_H
