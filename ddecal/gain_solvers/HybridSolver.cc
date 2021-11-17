// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "HybridSolver.h"

#include <cassert>

namespace dp3 {
namespace ddecal {

void HybridSolver::AddSolver(std::unique_ptr<SolverBase> solver) {
  if (!solvers_.empty() && (solver->NSolutionPolarizations() !=
                            solvers_.front().first->NSolutionPolarizations())) {
    throw std::runtime_error(
        "Solvers with different nr of polarizations can't be combined in the "
        "hybrid solver");
  }
  const size_t max_iterations = solver->GetMaxIterations();
  solvers_.emplace_back(std::move(solver), max_iterations);
}

HybridSolver::SolveResult HybridSolver::Solve(
    const SolveData& solve_data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(!solvers_.empty());
  size_t available_iterations = SolverBase::GetMaxIterations();
  SolveResult result;
  bool is_converged = false;
  for (const std::pair<std::unique_ptr<SolverBase>, size_t>& solver_info :
       solvers_) {
    solver_info.first->SetMaxIterations(solver_info.second);
    is_converged = RunSolver(*solver_info.first, available_iterations, result,
                             solve_data, solutions, time, stat_stream);
    if (is_converged && StopOnConvergence()) return result;
  }
  if (!is_converged) result.iterations = SolverBase::GetMaxIterations() + 1;
  return result;
}

bool HybridSolver::RunSolver(SolverBase& solver, size_t& available_iterations,
                             SolveResult& result, const SolveData& solve_data,
                             std::vector<std::vector<DComplex>>& solutions,
                             double time, std::ostream* stat_stream) {
  if (solver.GetMaxIterations() > available_iterations)
    solver.SetMaxIterations(available_iterations);
  SolveResult nextResult =
      solver.Solve(solve_data, solutions, time, stat_stream);
  result.iterations += nextResult.iterations;
  result.constraint_iterations += nextResult.constraint_iterations;
  result.results = std::move(nextResult.results);
  const bool is_converged = nextResult.iterations <= solver.GetMaxIterations();
  if (is_converged)
    available_iterations -= result.iterations;
  else  // If not converged, the solver has taken the maximum nr of
        // iterations.
    available_iterations -= solver.GetMaxIterations();
  return is_converged;
}

}  // namespace ddecal
}  // namespace dp3
