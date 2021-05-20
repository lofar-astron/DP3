// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "HybridSolver.h"

namespace dp3 {
namespace ddecal {

void HybridSolver::Initialize(size_t nAntennas, size_t nDirections,
                              size_t nChannels, size_t nChannelBlocks,
                              const std::vector<int>& ant1,
                              const std::vector<int>& ant2) {
  RegularSolverBase::Initialize(nAntennas, nDirections, nChannels,
                                nChannelBlocks, ant1, ant2);
  for (const std::pair<std::unique_ptr<RegularSolverBase>, size_t>&
           solver_info : SolversAndIterations()) {
    solver_info.first->Initialize(nAntennas, nDirections, nChannels,
                                  nChannelBlocks, ant1, ant2);
  }
}

}  // namespace ddecal
}  // namespace dp3
