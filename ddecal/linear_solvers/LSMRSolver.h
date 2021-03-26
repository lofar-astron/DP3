// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_LSMR_SOLVER_H
#define DDECAL_LSMR_SOLVER_H

#include <complex>
#include <vector>
#include <cmath>

#include "LLSSolver.h"

using complex = std::complex<float>;

namespace dp3 {
namespace base {

/* zLSMR prototype */
extern "C" void __clsmrmodule_MOD_clsmr(
    int* m, int* n, void Aprod1(int* m, int* n, complex* x, complex* y),
    void Aprod2(int* m, int* n, complex* x, complex* y), complex* b,
    float* damp, float* atol, float* btol, float* conlim, int* itnlim,
    int* localSize, int* nout, complex* x, int* istop, int* itn, float* normA,
    float* condA, float* normr, float* normAr, float* normx);

class LSMRSolver final : public LLSSolver {
 public:
  LSMRSolver(int m, int n, int nrhs);
  /**
   * Find X that minimizes || B - A*X || (i.e. best solution to A*X = B ; A
   * (MxN) * X (N) = B (M) Inputs are ordered column-major.
   * @param a The M-by-N matrix A
   * @param b On entry: input vector of size M, but with
   * dimension max(M, N) On succesful exit: the solution vector x.
   */

  bool Solve(complex* a, complex* b, double atolerance, double btolerance);

  virtual bool Solve(complex* a, complex* b) override;

  virtual bool Solve(complex* a, complex* b,
                     complex* x) override;  // use initial guess

  virtual void SetTolerance(const double tolerance) override {
    tolerance_ = tolerance;
  }

 private:
  static void Aprod1(int* m, int* n, complex* x, complex* y);
  static void Aprod2(int* m, int* n, complex* x, complex* y);

  thread_local static complex* LSMR_Matrix_A;

  double tolerance_;

  std::vector<complex> x_;
};

}  // namespace base
}  // namespace dp3

#endif
