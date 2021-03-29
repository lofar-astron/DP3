// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_ITERATIVE_SCALAR_SOLVER_H
#define DDECAL_ITERATIVE_SCALAR_SOLVER_H

#include "SolverBase.h"

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

/**
 * Implementation of a direction-dependent antenna scalar gain solver. It works
 * by iterating and assuming that, inside one iteration, antennas and directions
 * are independent. This results in a new solution for each antenna and
 * direction in each iteration, and after each iteration a small step is made
 * towards the gain.
 *
 * @ref IterativeDiagonalSolver implements the same algorithm, but solves for
 * diagonal gains.
 *
 * It is a counterpart of @ref ScalarSolver, but faster, possibly at the expense
 * of stability, but this has not been properly compared.
 */
class IterativeScalarSolver final : public SolverBase {
 public:
  IterativeScalarSolver() : SolverBase() {}

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
  void AddOrSubtractDirection(std::vector<std::vector<Complex>>& v_residual,
                              size_t channel_block_index, size_t direction,
                              const std::vector<DComplex>& solutions);

  void SolveDirection(size_t channel_block_index,
                      const std::vector<std::vector<Complex>>& v_residual,
                      size_t direction, const std::vector<DComplex>& solutions,
                      std::vector<DComplex>& next_solutions);
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_ITERATIVE_SCALAR_SOLVER_H
