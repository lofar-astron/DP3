// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LSTSQ_SOLVER_H
#define LSTSQ_SOLVER_H

#include <boost/algorithm/string/case_conv.hpp>

#include <cmath>
#include <complex>
#include <vector>
#include <memory>

namespace DP3 {
namespace DPPP {

/**
 * Linear least-squares solver to use.
 */
enum class LLSSolverType { QR, SVD, LSMR };

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
    if (boost::algorithm::to_lower_copy(solver) == "lsmr") {
      return LLSSolverType::LSMR;
    } else if (boost::algorithm::to_lower_copy(solver) == "svd") {
      return LLSSolverType::SVD;
    } else if (boost::algorithm::to_lower_copy(solver) == "qr") {
      return LLSSolverType::QR;
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

}  // namespace DPPP
}  // namespace DP3

#endif
