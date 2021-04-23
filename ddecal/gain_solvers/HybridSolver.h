// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_HYBRID_SOLVER_H
#define DDECAL_HYBRID_SOLVER_H

#include "SolverBase.h"
#include "DiagonalSolver.h"
#include "IterativeDiagonalSolver.h"

namespace dp3 {
namespace base {

class HybridSolver final : public SolverBase {
 public:
  HybridSolver() : SolverBase(), stop_on_convergence_(false) {}

  SolveResult Solve(const SolverBuffer& solver_buffer,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  void Initialize(size_t nAntennas, size_t nDirections, size_t nChannels,
                  size_t nChannelBlocks, const std::vector<int>& ant1,
                  const std::vector<int>& ant2) override;

  size_t NSolutionPolarizations() const override {
    return solvers_.empty() ? 0
                            : solvers_.front().first->NSolutionPolarizations();
  }

  void AddSolver(std::unique_ptr<SolverBase> solver);

  void SetStopOnConvergence(bool stop_on_convergence) {
    stop_on_convergence_ = stop_on_convergence;
  }

 private:
  // List of solvers with their maximum number of iterations
  std::vector<std::pair<std::unique_ptr<SolverBase>, size_t>> solvers_;
  bool stop_on_convergence_;

  static bool RunSolver(SolverBase& solver, size_t& available_iterations,
                        SolveResult& result, const SolverBuffer& solver_buffer,
                        std::vector<std::vector<DComplex>>& solutions,
                        double time, std::ostream* stat_stream);
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_HYBRID_SOLVER_H
