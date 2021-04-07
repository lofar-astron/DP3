// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_FULL_JONES_SOLVER_H
#define DDECAL_FULL_JONES_SOLVER_H

#include "SolverBase.h"

#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class FullJonesSolver final : public SolverBase {
 public:
  FullJonesSolver() : SolverBase() {}

  /**
   * Solves full Jones matrices.
   * @param data_buffers as in @ref SolverBase::Solve()
   * @param model_buffers as in @ref SolverBase::Solve()
   * @param solutions An array, where @c solutions[ch] is a pointer to
   * channelblock @c ch, that points to antenna x directions solutions. Each
   * solution consists of 4 complex values forming the full Jones matrix.
   */
  SolveResult Solve(const std::vector<DPBuffer>& unweighted_data_buffers,
                    const std::vector<std::vector<DPBuffer*>>& model_buffers,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

 private:
  void PerformIteration(
      const std::vector<std::vector<DPBuffer*>>& model_buffers,
      size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
      std::vector<Matrix>& vs, const std::vector<DComplex>& solutions,
      std::vector<DComplex>& nextSolutions);
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_FULL_JONES_SOLVER_H
