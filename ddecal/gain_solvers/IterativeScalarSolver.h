// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_ITERATIVE_SCALAR_SOLVER_H
#define DDECAL_ITERATIVE_SCALAR_SOLVER_H

#include "RegularSolverBase.h"

#include <complex>
#include <vector>

namespace dp3 {
namespace ddecal {

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
class IterativeScalarSolver final : public RegularSolverBase {
 public:
  IterativeScalarSolver() : RegularSolverBase() {}

  SolveResult Solve(const SolverBuffer& solver_buffer,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 1; }

 private:
  void PerformIteration(const SolverBuffer& solver_buffer,
                        size_t channelBlockIndex,
                        std::vector<std::vector<Complex>>& v_residual,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions);

  template <bool Add>
  void AddOrSubtractDirection(const SolverBuffer& solver_buffer,
                              std::vector<std::vector<Complex>>& v_residual,
                              size_t channel_block_index, size_t direction,
                              const std::vector<DComplex>& solutions);

  void SolveDirection(const SolverBuffer& solver_buffer,
                      size_t channel_block_index,
                      const std::vector<std::vector<Complex>>& v_residual,
                      size_t direction, const std::vector<DComplex>& solutions,
                      std::vector<DComplex>& next_solutions);
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_ITERATIVE_SCALAR_SOLVER_H
