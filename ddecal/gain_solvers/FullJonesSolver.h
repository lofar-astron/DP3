// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_FULL_JONES_SOLVER_H
#define DDECAL_FULL_JONES_SOLVER_H

#include "SolverBase.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

class FullJonesSolver final : public SolverBase {
 public:
  SolveResult Solve(const SolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 4; }

 private:
  void PerformIteration(const SolveData::ChannelBlockData& cb_data,
                        std::vector<Matrix>& g_times_cs,
                        std::vector<Matrix>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& next_solutions);
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_BDA_FULL_JONES_SOLVER_H
