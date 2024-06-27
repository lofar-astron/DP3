// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_ITERATIVE_SCALAR_SOLVER_H
#define DDECAL_ITERATIVE_SCALAR_SOLVER_H

#include "SolverBase.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

/**
 * Iterative scalar solver supporting direction-dependent solution intervals.
 */
class IterativeScalarSolver final : public SolverBase {
 public:
  SolveResult Solve(const FullSolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 1; }

  bool SupportsDdSolutionIntervals() const override { return true; }

 private:
  void PerformIteration(size_t ch_block,
                        const FullSolveData::ChannelBlockData& cb_data,
                        std::vector<aocommon::MC2x2F>& v_residual,
                        const std::vector<DComplex>& solutions,
                        SolutionTensor& next_solutions);

  template <bool Add>
  void AddOrSubtractDirection(const FullSolveData::ChannelBlockData& cb_data,
                              std::vector<aocommon::MC2x2F>& v_residual,
                              size_t direction,
                              const std::vector<DComplex>& solutions);

  void SolveDirection(size_t ch_block,
                      const FullSolveData::ChannelBlockData& cb_data,
                      const std::vector<aocommon::MC2x2F>& v_residual,
                      size_t direction, const std::vector<DComplex>& solutions,
                      SolutionTensor& next_solutions);
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_ITERATIVE_SCALAR_SOLVER_H
