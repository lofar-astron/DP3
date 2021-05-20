// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_HYBRID_SOLVER_H
#define DDECAL_BDA_HYBRID_SOLVER_H

#include "HybridSolverBase.h"
#include "BdaSolverBase.h"

namespace dp3 {
namespace ddecal {

class BdaHybridSolver final : public HybridSolverBase<BdaSolverBase> {
 public:
  SolveResult Solve(const SolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override {
    return CallSolvers(data, solutions, time, stat_stream);
  }

  void Initialize(size_t n_antennas, size_t n_directions,
                  size_t n_channel_blocks) override {
    BdaSolverBase::Initialize(n_antennas, n_directions, n_channel_blocks);
    for (const std::pair<std::unique_ptr<BdaSolverBase>, size_t>& solver_info :
         SolversAndIterations()) {
      solver_info.first->Initialize(n_antennas, n_directions, n_channel_blocks);
    }
  }
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_HYBRID_SOLVER_H
