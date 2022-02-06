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
              const SolverMode mode)
      : robust_nu_(robust_nu),
        batch_iter_(max_iter),
        history_size_(history_size),
        minibatches_(minibatches),
        mode_(mode){};

  SolveResult Solve(const SolveData& data,
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
  SolverMode GetMode() const { return mode_; }

 private:
  double robust_nu_;
  size_t batch_iter_;
  size_t history_size_;
  size_t minibatches_;
  SolverMode mode_;
};
#endif /* HAVE_LIBDIRAC */
}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_LBFGS_SOLVER_H
