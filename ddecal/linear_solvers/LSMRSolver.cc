// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef USE_LSMR

#include "LSMRSolver.h"

#include <algorithm>
#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

namespace dp3 {
namespace ddecal {

thread_local complex* LSMRSolver::LSMR_Matrix_A;

LSMRSolver::LSMRSolver(int m, int n, int nrhs)
    : LLSSolver(m, n, nrhs), tolerance_(1.0E-7), x_(n) {}

void LSMRSolver::Aprod1(int* m, int* n, complex* x, complex* y) {
  // compute y:= y + A x
  complex* A = LSMR_Matrix_A;
  for (int mp = 0; mp < *m; ++mp) {
    for (int np = 0; np < *n; ++np) {
      y[mp] += A[mp + np * *m] * x[np];  // A is stored in column-major order
    }
  }
}
void LSMRSolver::Aprod2(int* m, int* n, complex* x, complex* y) {
  // compute x := x + A^T y
  complex* A = LSMR_Matrix_A;
  for (int np = 0; np < *n; ++np) {
    for (int mp = 0; mp < *m; ++mp) {
      x[np] += std::conj(A[np * *m + mp]) * y[mp];
    }
  }
}

/**
 * Find X that minimizes || B - A*X || (i.e. best solution to A*X = B ; A
 * (MxN) * X (N) = B (M) Inputs are ordered column-major.
 * @param a The M-by-N matrix A
 * @param b On entry: input vector of size M, but with
 * dimension max(M, N) On succesful exit: the solution vector x.
 */

bool LSMRSolver::Solve(complex* a, complex* b, double atolerance,
                       double btolerance) {
  LSMR_Matrix_A = a;

  float damp = 0.0;
  float atol = atolerance;
  float btol = btolerance;
  float conlim = 1.0E8;
  int itnlim = 4 * n_;
  int localSize = std::min(n_, m_);
  int nout = 0;
  int istop;
  int itn;
  float normA;
  float condA;
  float normr;
  float normAr;
  float normx;

  __clsmrmodule_MOD_clsmr(&m_, &n_, Aprod1, Aprod2, b, &damp, &atol, &btol,
                          &conlim, &itnlim, &localSize, &nout, x_.data(),
                          &istop, &itn, &normA, &condA, &normr, &normAr,
                          &normx);

  std::copy_n(x_.data(), n_, b);

  //  for (int i = 0; i < n_; ++i) {
  //    b[i] = x_[i];
  // }

  return itn < itnlim;
}

bool LSMRSolver::Solve(complex* a, complex* b) {
  bool success = Solve(a, b, tolerance_, tolerance_);
  return success;
}

bool LSMRSolver::Solve(complex* a, complex* b, complex* x) {
  double normb = 0.0;
  for (int i = 0; i < m_; ++i) {
    normb += (double)(b[i].real() * b[i].real() + b[i].imag() * b[i].imag());
  }
  normb = std::sqrt(normb);

  for (int np = 0; np < n_; ++np) {
    x[np] *= -1.0;  // -x
  }

  LSMR_Matrix_A = a;
  Aprod1(&m_, &n_, x, b);  // b - Ax

  // compute norm of r0 to use in calculating tolerance
  double normr0 = 0.0;
  for (int i = 0; i < m_; ++i) {
    normr0 += (double)(b[i].real() * b[i].real() + b[i].imag() * b[i].imag());
  }
  normr0 = std::sqrt(normr0);

  double atolerance = tolerance_;
  double btolerance = tolerance_ * normb / normr0;

  bool success = Solve(a, b, atolerance, btolerance);

  for (int np = 0; np < n_; ++np) {
    b[np] -= x[np];  // result + x0, x was inverted above
  }

  return success;
}

}  // namespace ddecal
}  // namespace dp3

#endif
