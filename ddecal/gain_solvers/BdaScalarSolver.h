// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_SCALAR_SOLVER_H
#define DDECAL_BDA_SCALAR_SOLVER_H

#include "BdaSolverBase.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

class BdaScalarSolver final : public BdaSolverBase {
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
                        std::vector<DComplex>& next_solutions,
                        double iteration_fraction, double solver_precision);
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_SCALAR_SOLVER_H
