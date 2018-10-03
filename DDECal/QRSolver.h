#ifndef QR_SOLVER_H
#define QR_SOLVER_H

#include <complex>
#include <vector>
#include <cmath>

typedef std::complex<double> dcomplex;

/* ZGELS prototype */
extern "C" void zgels_(char* trans, int* m, int* n, int* nrhs, dcomplex* a,
  int* lda, dcomplex* b, int* ldb, dcomplex* work, int* lwork, int* info);

class QRSolver
{
public:
  QRSolver(int m, int n, int nrhs) :
    _m(m), _n(n), _nrhs(nrhs)
  { }
  
  /**
   * Find X that minimizes || B - A*X || (i.e. best solution to A*X = B ; A (MxN) * X (NxNRHS) = B (MxNRHS)
   * Inputs are ordered column-major.
   * @param a The M-by-N matrix A
   * @param b On entry: input matrix of size M x NRHS, but with leading dimension max(M, N)
   * On succesful exit: the solution vectors, stored columnwise (N, NRHS)
   */
  bool Solve(dcomplex* a, dcomplex* b)
  {
    int info;
    char trans = 'N'; // No transpose
    int ldb = std::max(_m, _n);
    if(_work.empty())
    {
      int lwork = -1;
      dcomplex wkopt;
      zgels_(&trans, &_m, &_n, &_nrhs, a, &_m, b, &ldb, &wkopt, &lwork, &info );
      _work.resize((int) wkopt.real());
    }
    int lwork = _work.size();
    zgels_(&trans, &_m, &_n, &_nrhs, a, &_m, b, &ldb, _work.data(), &lwork, &info);
    // Check for full rank
    return info == 0;
  }
private:
  /** Number of rows in matrix A */
  int _m;
  
  /** Number of cols in matrix A */
  int _n;
  
  /** Number of columns in matrices B and X */
  int _nrhs;
  
  std::vector<dcomplex> _work;
};

#endif
