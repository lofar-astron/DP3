// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_ITERATIVE_DIAGONAL_SOLVER_H
#define DDECAL_ITERATIVE_DIAGONAL_SOLVER_H

#include "SolverBase.h"

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class IterativeDiagonalSolver final : public SolverBase {
 public:
  IterativeDiagonalSolver() : SolverBase() {}

  SolveResult Solve(const std::vector<Complex*>& unweighted_data,
                    const std::vector<float*>& weights,
                    std::vector<std::vector<Complex*>>&& unweighted_model_data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

 private:
  void PerformIteration(size_t channelBlockIndex,
                        std::vector<std::vector<Complex>>& v_residual,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions);

  template <bool Add>
  void SubtractDirection(std::vector<std::vector<Complex>>& v_residual,
                         size_t channel_block_index, size_t direction,
                         const std::vector<DComplex>& solutions);

  void SolveDirection(size_t channel_block_index,
                      const std::vector<std::vector<Complex>>& v_residual,
                      size_t direction, const std::vector<DComplex>& solutions,
                      std::vector<DComplex>& next_solutions);
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_ITERATIVE_DIAGONAL_SOLVER_H
