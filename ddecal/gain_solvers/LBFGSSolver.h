// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_LBFGS_SOLVER_H
#define DDECAL_LBFGS_SOLVER_H

#include "SolverBase.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

#ifdef HAVE_LIBDIRAC
class LBFGSSolver final : public SolverBase {
 public:
  enum SolverMode { kFull, kDiagonal, kScalar };
  LBFGSSolver(const double robust_nu, const size_t max_iter,
              const size_t history_size, const size_t minibatches,
              const double min_solution, const double max_solution,
              const SolverMode mode)
      : robust_nu_(robust_nu),
        batch_iter_(max_iter),
        history_size_(history_size),
        minibatches_(minibatches),
        min_solution_(min_solution),
        max_solution_(max_solution),
        bound_constrained_((min_solution || max_solution) &&
                           (min_solution != max_solution)),
        mode_(mode){};

  /// Split real and imaginary parts and combine them into a single tensor.
  /// This internal function is static, since it has unit tests.
  /// @param solutions A vector with complex double values.
  /// @return An XTensor with double values, which has twice the length of
  ///         the solutions argument. The first and second halves of the
  ///         tensor will contain the real and imaginary parts of the solutions,
  ///         respectively.
  static xt::xtensor<double, 1> SplitSolutions(
      const std::vector<DComplex>& solutions);

  /// Merge the real and imaginary part of d_storage into next_solutions.
  /// This internal function is static, since it has unit tests.
  /// @param next_solutions [out] Target for writing the solutions.
  /// @param ch_block Channel block index for next_solutions.
  /// @param d_storage An 1D xtensor. The first and second halves contain the
  ///                  real and imaginary parts of the solutions, respectively.
  static void MergeSolutions(SolutionTensor& next_solutions, size_t ch_block,
                             const xt::xtensor<double, 1>& d_storage);

  SolveResult Solve(const FullSolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override {
    switch (mode_) {
      case kFull:
        return 4;
      case kDiagonal:
        return 2;
      case kScalar:
        return 1;
    }
    assert(false);
    return 0;  // will not get here
  }

  double GetRobustDOF() const { return robust_nu_; }
  size_t GetMaxIter() const { return batch_iter_; }
  size_t GetHistorySize() const { return history_size_; }
  size_t GetMinibatches() const { return minibatches_; }
  double GetMinSolution() const { return min_solution_; }
  double GetMaxSolution() const { return max_solution_; }
  bool GetBoundConstrained() const { return bound_constrained_; }
  SolverMode GetMode() const { return mode_; }

 private:
  double robust_nu_;
  size_t batch_iter_;
  size_t history_size_;
  size_t minibatches_;
  double min_solution_;
  double max_solution_;
  bool bound_constrained_;
  SolverMode mode_;
};
#endif /* HAVE_LIBDIRAC */
}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_LBFGS_SOLVER_H
