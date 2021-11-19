// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_LLS_SOLVER_H
#define DDECAL_LLS_SOLVER_H

#include <complex>
#include <memory>
#include <utility>
#include <vector>

namespace dp3 {
namespace ddecal {

/**
 * Linear least-squares solver to use.
 */
enum class LLSSolverType { QR, SVD, NORMAL_EQUATIONS };

/**
 * Base class for linear least-square solvers. These are used by
 * the solution solvers such as ScalarSolver, DiagonalSolver and
 * FullJonesSolver.
 */
class LLSSolver {
 public:
  using complex = std::complex<float>;
  LLSSolver(int m, int n, int nrhs) : m_(m), n_(n), nrhs_(nrhs) {}

  virtual ~LLSSolver() {}

  virtual bool Solve(complex* a, complex* b) = 0;
  virtual bool Solve(complex* a, complex* b,
                     [[maybe_unused]] complex* initial_value) {
    return Solve(a, b);
  }  // can use initial guess x

  static std::unique_ptr<LLSSolver> Make(LLSSolverType lss_type, int m, int n,
                                         int nrhs);

  static LLSSolverType ParseType(const std::string& solver);

 protected:
  int m_;     // number of rows in matrix A
  int n_;     // number of columns in matrix A
  int nrhs_;  // number of columns in vector b
};

}  // namespace ddecal
}  // namespace dp3

#endif
