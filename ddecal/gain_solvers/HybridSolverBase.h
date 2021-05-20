// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_HYBRID_SOLVER_BASE_H
#define DDECAL_HYBRID_SOLVER_BASE_H

#include "SolverBase.h"

#include <cassert>

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
 * @ref HybridSolver and @ref BdaHybridSolver make use of this class.
 *
 * @tparam Base Base solver class for the internal solvers. It should be either
 * RegularSolver or BdaSolver. This template argument allows building a generic
 * hybrid solver that supports both regular and BDA solvers.
 */
template <typename Base>
class HybridSolverBase : public Base {
 public:
  HybridSolverBase() : Base(), stop_on_convergence_(false) {}

  size_t NSolutionPolarizations() const override {
    return solvers_.empty() ? 0
                            : solvers_.front().first->NSolutionPolarizations();
  }

  void SetNThreads(size_t n_threads) override {
    Base::SetNThreads(n_threads);
    for (std::pair<std::unique_ptr<Base>, size_t>& solver_info : solvers_) {
      solver_info.first->SetNThreads(n_threads);
    }
  }

  /**
   * Add a solver. Solvers should be added in the order that they should
   * be called, and each solver should have its maximum number of iterations set
   * before adding it.
   */
  void AddSolver(std::unique_ptr<Base> solver) {
    if (!solvers_.empty()) {
      if (solver->NSolutionPolarizations() !=
          solvers_.front().first->NSolutionPolarizations())
        throw std::runtime_error(
            "Solvers with different nr of polarizations can't be combined in "
            "the hybrid solver");
    }
    size_t iter = solver->GetMaxIterations();
    solvers_.emplace_back(std::move(solver), iter);
  }

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

 protected:
  using typename Base::DComplex;
  using typename Base::SolveResult;

  const std::vector<std::pair<std::unique_ptr<Base>, size_t>>&
  SolversAndIterations() {
    return solvers_;
  }

  template <typename SolveBufferType>
  SolveResult CallSolvers(const SolveBufferType& solver_buffer,
                          std::vector<std::vector<DComplex>>& solutions,
                          double time, std::ostream* stat_stream) {
    assert(!solvers_.empty());
    size_t available_iterations = Base::GetMaxIterations();
    SolveResult result;
    bool is_converged = false;
    for (const std::pair<std::unique_ptr<Base>, size_t>& solver_info :
         solvers_) {
      solver_info.first->SetMaxIterations(solver_info.second);
      is_converged = RunSolver(*solver_info.first, available_iterations, result,
                               solver_buffer, solutions, time, stat_stream);
      if (is_converged && StopOnConvergence()) return result;
    }
    if (!is_converged) result.iterations = Base::GetMaxIterations() + 1;
    return result;
  }

 private:
  template <typename SolveBufferType>
  bool RunSolver(Base& solver, size_t& available_iterations,
                 SolveResult& result, const SolveBufferType& solver_buffer,
                 std::vector<std::vector<DComplex>>& solutions, double time,
                 std::ostream* stat_stream) {
    if (solver.GetMaxIterations() > available_iterations)
      solver.SetMaxIterations(available_iterations);
    SolveResult nextResult =
        solver.Solve(solver_buffer, solutions, time, stat_stream);
    result.iterations += nextResult.iterations;
    result.constraint_iterations += nextResult.constraint_iterations;
    result.results = std::move(nextResult.results);
    const bool is_converged =
        nextResult.iterations <= solver.GetMaxIterations();
    if (is_converged)
      available_iterations -= result.iterations;
    else  // If not converged, the solver has taken the maximum nr of
          // iterations.
      available_iterations -= solver.GetMaxIterations();
    return is_converged;
  }

  // List of solvers with their maximum number of iterations
  std::vector<std::pair<std::unique_ptr<Base>, size_t>> solvers_;
  bool stop_on_convergence_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
