// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef NORMALEQ_SOLVER_H
#define NORMALEQ_SOLVER_H

#include "LLSSolver.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace dp3 {
namespace ddecal {

/* CGELSD prototype */
extern "C" void cposv_(char* uplo, int* n, int* nrhs, std::complex<float>* a,
                       int* lda, std::complex<float>* b, int* ldb, int* info);

class NormalEquationsSolver final : public LLSSolver {
 public:
  NormalEquationsSolver(int m, int n, int nrhs)
      : LLSSolver(m, n, nrhs), adaggera_(n_ * n_), adaggerb_(n_ * nrhs_) {}

  /**
   * Find X that minimizes || B - A*X || (i.e. best solution to A*X = B ; A
   * (MxN) * X (NxNRHS) = B (MxNRHS) Inputs are ordered column-major.
   * @param a The M-by-N matrix A
   * @param b On entry: input matrix of size M x NRHS, but with leading
   * dimension max(M, N) On succesful exit: the solution vectors, stored
   * columnwise (N, NRHS)
   */
  bool Solve(std::complex<float>* a, std::complex<float>* b) override {
    // compute a^dagger a
    for (int n = 0; n < n_; ++n) {
      for (int np = n; np < n_; ++np) {  // fill upper triangle only
        adaggera_[n + np * n_] = std::complex<float>(0.0, 0.0);
        for (int m = 0; m < m_; ++m) {
          adaggera_[n + np * n_] += std::conj(a[m + n * m_]) * a[m + np * m_];
        }
      }
    }
    for (int p = 0; p < nrhs_; ++p) {
      for (int n = 0; n < n_; ++n) {
        adaggerb_[n + p * n_] = std::complex<float>(0.0, 0.0);
        for (int m = 0; m < m_; ++m) {
          adaggerb_[n + p * n_] += std::conj(a[m + n * m_]) * b[m + p * m_];
        }
      }
    }

    int info;
    char uplo = 'U';
    int ldb = n_;

    // solve Hermitian system of normal equations using Cholesky decomposition
    cposv_(&uplo, &n_, &nrhs_, adaggera_.data(), &n_, adaggerb_.data(), &ldb,
           &info);

    std::copy_n(adaggerb_.data(), n_ * nrhs_, b);

    // Check for full rank
    return info == 0;
  }

 private:
  std::vector<std::complex<float>> adaggera_;
  std::vector<std::complex<float>> adaggerb_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
