// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef QR_SOLVER_H
#define QR_SOLVER_H

#include <complex>
#include <vector>
#include <cmath>

#include "LLSSolver.h"

namespace dp3 {
namespace ddecal {

using complex = std::complex<float>;

/* CGELS prototype */
extern "C" void cgels_(char* trans, int* m, int* n, int* nrhs, complex* a,
                       int* lda, complex* b, int* ldb, complex* work,
                       int* lwork, int* info);

class QRSolver final : public LLSSolver {
 public:
  QRSolver(int m, int n, int nrhs) : LLSSolver(m, n, nrhs) {}

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
    char trans = 'N';  // No transpose
    int ldb = std::max(m_, n_);
    if (_work.empty()) {
      int lwork = -1;
      complex wkopt;
      cgels_(&trans, &m_, &n_, &nrhs_, a, &m_, b, &ldb, &wkopt, &lwork, &info);
      _work.resize((int)wkopt.real());
    }
    int lwork = _work.size();
    cgels_(&trans, &m_, &n_, &nrhs_, a, &m_, b, &ldb, _work.data(), &lwork,
           &info);
    /// Check for full rank
    return info == 0;
  }

 private:
  std::vector<complex> _work;
};

}  // namespace ddecal
}  // namespace dp3

#endif
