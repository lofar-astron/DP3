// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_HYBRID_SOLVER_H
#define DDECAL_HYBRID_SOLVER_H

#include "SolverBase.h"

namespace dp3 {
namespace ddecal {

/**
 * This class implements the functionality for a hybrid solver. A hybrid solver
 * maintains a list of solvers of different types and calls them one by one
 * with a provided number of maximum iterations.
 *
 * All solvers must solve the same number of solutions. The solvers can have
 * different settings and constraints.
 *
 * This allows starting with a slow, stable solver to converge when the initial
 * values are quite off, and switches to faster solvers when approaching the
 * correct solutions.
 */
class HybridSolver final : public SolverBase {
 public:
  HybridSolver() : stop_on_convergence_(false) {}

  void Initialize(size_t n_antennas,
                  const std::vector<size_t>& n_solutions_per_direction,
                  size_t n_channel_blocks) override {
    SolverBase::Initialize(n_antennas, n_solutions_per_direction,
                           n_channel_blocks);
    for (const std::pair<std::unique_ptr<SolverBase>, size_t>& solver_info :
         solvers_) {
      solver_info.first->Initialize(n_antennas, n_solutions_per_direction,
                                    n_channel_blocks);
    }
  }

  size_t NSolutionPolarizations() const override {
    return solvers_.empty() ? 0
                            : solvers_.front().first->NSolutionPolarizations();
  }

  void SetNThreads(size_t n_threads) override {
    SolverBase::SetNThreads(n_threads);
    for (std::pair<std::unique_ptr<SolverBase>, size_t>& solver_info :
         solvers_) {
      solver_info.first->SetNThreads(n_threads);
    }
  }

  /**
   * Add a solver. Solvers should be added in the order that they should
   * be called, and each solver should have its maximum number of iterations set
   * before adding it.
   */
  void AddSolver(std::unique_ptr<SolverBase> solver);

  /**
   * List of solvers that need constraint initialization. This
   * list does not include the Hybrid(Base) solver, i.e., it does
   * not contain the "this" pointer.
   */
  std::vector<SolverBase*> ConstraintSolvers() override {
    std::vector<SolverBase*> solvers;
    for (const auto& solverinfo : solvers_)
      solvers.push_back(solverinfo.first.get());
    return solvers;
  }

  /**
   * @{
   * If set to true and a solver finishes before reaching its maximum number
   * of iterations, subsequent solvers are not called.
   */
  bool StopOnConvergence() const { return stop_on_convergence_; }
  void SetStopOnConvergence(bool stop_on_convergence) {
    stop_on_convergence_ = stop_on_convergence;
  }
  /** @} */

  SolveResult Solve(const SolveData& solve_data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

 private:
  bool RunSolver(SolverBase& solver, size_t& available_iterations,
                 SolveResult& result, const SolveData& solve_data,
                 std::vector<std::vector<DComplex>>& solutions, double time,
                 std::ostream* stat_stream);

  // List of solvers with their maximum number of iterations
  std::vector<std::pair<std::unique_ptr<SolverBase>, size_t>> solvers_;
  bool stop_on_convergence_;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_HYBRID_SOLVER_H
