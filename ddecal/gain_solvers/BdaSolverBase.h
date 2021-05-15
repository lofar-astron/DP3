// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDE_BDA_SOLVER_BASE_H
#define DDE_BDA_SOLVER_BASE_H

#include "SolverBase.h"

#include "../constraints/Constraint.h"

#include "../linear_solvers/LLSSolver.h"

#include <boost/algorithm/string.hpp>

#include <cassert>
#include <complex>
#include <vector>
#include <memory>

namespace dp3 {
namespace ddecal {

class SolveData;

class BdaSolverBase : public SolverBase {
 public:
  /**
   * Solves multi-directional Jones matrices. Takes the (single) measured data
   * and the (multi-directional) model data, and solves the optimization
   * problem that minimizes the norm of the differences.
   *
   * @param data Buffer with weighted data and model data.
   * @param solutions The per-channel and per-antenna solutions.
   * solutions[ch] is a pointer for channelblock ch to antenna x directions x
   * pol solutions.
   * @param statStream Optional pointer to a stream for displaying statistics.
   */
  virtual SolveResult Solve(const SolveData& data,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) = 0;
};

}  // namespace ddecal
}  // namespace dp3

#endif
