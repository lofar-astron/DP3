// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SCALAR_SOLVER_H
#define DDECAL_SCALAR_SOLVER_H

#include "SolverBase.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

class ScalarSolver final : public SolverBase {
 public:
  SolveResult Solve(const SolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 1; }

 private:
  void PerformIteration(const SolveData::ChannelBlockData& cb_data,
                        std::vector<Matrix>& g_times_cs,
                        std::vector<Matrix>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& next_solutions);

  /**
   * Initialize the model matrix for a channel block.
   *
   * The number of elements in the model matrix depends on the number of
   * antenna visibilities in the corresponding channel block.
   */
  void InitializeModelMatrix(
      const SolveData::ChannelBlockData& channel_block_data,
      std::vector<Matrix>& g_times_cs, std::vector<Matrix>& vs) const;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_SCALAR_SOLVER_H
