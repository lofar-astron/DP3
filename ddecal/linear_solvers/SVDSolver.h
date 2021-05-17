// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SVD_SOLVER_H
#define SVD_SOLVER_H

#include <complex>
#include <vector>
#include <cmath>

#include "LLSSolver.h"

namespace dp3 {
namespace ddecal {

using complex = std::complex<float>;

/* CGELSD prototype */
extern "C" void cgelss_(int* m, int* n, int* nrhs, complex* a, int* lda,
                        complex* b, int* ldb, float* s, float* rcond, int* rank,
                        complex* work, int* lwork, float* rwork, int* info);

class SVDSolver final : public LLSSolver {
 public:
  SVDSolver(int m, int n, int nrhs) : LLSSolver(m, n, nrhs) {}

  /**
   * Find X that minimizes || B - A*X || (i.e. best solution to A*X = B ; A
   * (MxN) * X (NxNRHS) = B (MxNRHS) Inputs are ordered column-major.
   * @param a The M-by-N matrix A
   * @param b On entry: input matrix of size M x NRHS, but with leading
   * dimension max(M, N) On succesful exit: the solution vectors, stored
   * columnwise (N, NRHS)
   */
  virtual bool Solve(complex* a, complex* b) override {
    int info;
    int ldb = std::max(m_, n_);
    std::vector<float> s;
    s.resize(std::min(n_, m_));
    float rcond = 0.0;
    int rank;
    std::vector<float> rwork(5 * std::min(m_, n_));

    if (work_.empty()) {
      int lwork = -1;
      complex wkopt;
      cgelss_(&m_, &n_, &nrhs_, a, &m_, b, &ldb, s.data(), &rcond, &rank,
              &wkopt, &lwork, rwork.data(), &info);
      work_.resize((int)wkopt.real());
    }
    int lwork = work_.size();
    cgelss_(&m_, &n_, &nrhs_, a, &m_, b, &ldb, s.data(), &rcond, &rank,
            work_.data(), &lwork, rwork.data(), &info);
    /// Check for full rank
    return info == 0;
  }

 private:
  std::vector<complex> work_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
