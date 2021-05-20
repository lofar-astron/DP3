// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_HYBRID_SOLVER_H
#define DDECAL_HYBRID_SOLVER_H

#include "HybridSolverBase.h"
#include "RegularSolverBase.h"

namespace dp3 {
namespace ddecal {

class HybridSolver final : public HybridSolverBase<RegularSolverBase> {
 public:
  SolveResult Solve(const SolverBuffer& solver_buffer,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override {
    return CallSolvers(solver_buffer, solutions, time, stat_stream);
  }

  void Initialize(size_t nAntennas, size_t nDirections, size_t nChannels,
                  size_t nChannelBlocks, const std::vector<int>& ant1,
                  const std::vector<int>& ant2) override;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_HYBRID_SOLVER_H
