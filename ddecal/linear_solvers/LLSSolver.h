// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_LLS_SOLVER_H
#define DDECAL_LLS_SOLVER_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <memory>

namespace dp3 {
namespace base {

/**
 * Linear least-squares solver to use.
 */
enum class LLSSolverType { QR, SVD, LSMR, NORMAL_EQUATIONS };

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

  virtual void SetTolerance([[maybe_unused]] double tolerance) {}

  static std::unique_ptr<LLSSolver> Make(LLSSolverType lss_type, int m, int n,
                                         int nrhs);

  static LLSSolverType ParseType(const std::string& solver) {
    std::string solver_lowercase = solver;
    std::transform(solver.begin(), solver.end(), solver_lowercase.begin(),
                   [](char c) { return std::tolower(c); });

    if (solver_lowercase == "lsmr") {
      return LLSSolverType::LSMR;
    } else if (solver_lowercase == "svd") {
      return LLSSolverType::SVD;
    } else if (solver_lowercase == "qr") {
      return LLSSolverType::QR;
    } else if (solver_lowercase == "normalequations") {
      return LLSSolverType::NORMAL_EQUATIONS;
    } else {
      throw std::runtime_error("Unknown least squares solver requested: " +
                               solver);
    }
  }

  static std::pair<double, double> ParseTolerances(double min_tolerance,
                                                   double max_tolerance) {
    double result_min = min_tolerance;
    double result_max = max_tolerance;
    if (result_max < result_min) {
      result_max = result_min;
    }
    return std::pair<double, double>(result_min, result_max);
  }

 protected:
  int m_;     // number of rows in matrix A
  int n_;     // number of columns in matrix A
  int nrhs_;  // number of columns in vector b
};

}  // namespace base
}  // namespace dp3

#endif
