// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_GAIN_SOLVERS_KERNELS_MATRIXCOMPLEX2X2_H_
#define DP3_DDECAL_GAIN_SOLVERS_KERNELS_MATRIXCOMPLEX2X2_H_

#include <cuComplex.h>

template <typename T>
struct cuM2x2 {
 public:
  __device__ cuM2x2() {
    data[0] = {0};
    data[1] = {0};
    data[2] = {0};
    data[3] = {0};
  }
  __device__ cuM2x2(T a, T b, T c, T d) {
    data[0] = a;
    data[1] = b;
    data[2] = c;
    data[3] = d;
  }

  inline __device__ const T& operator[](int i) const { return data[i]; }
  inline __device__ T& operator[](int i) { return data[i]; }

 private:
  T data[4];
};

template <typename T>
struct cuM2x2Diagonal {
 public:
  __device__ cuM2x2Diagonal() {
    data[0] = {0};
    data[1] = {0};
  }
  __device__ cuM2x2Diagonal(T a, T b) {
    data[0] = a;
    data[1] = b;
  }

  inline __device__ const T& operator[](int i) const { return data[i]; }
  inline __device__ T& operator[](int i) { return data[i]; }

 private:
  T data[2];
};

using cuM2x2FloatComplex = cuM2x2<cuFloatComplex>;
using cuM2x2DoubleComplex = cuM2x2<cuDoubleComplex>;

inline __device__ cuM2x2FloatComplex operator+(const cuM2x2FloatComplex& a,
                                               const cuM2x2FloatComplex& b) {
  return cuM2x2FloatComplex(cuCaddf(a[0], b[0]), cuCaddf(a[1], b[1]),
                            cuCaddf(a[2], b[2]), cuCaddf(a[3], b[3]));
}

inline __device__ cuM2x2FloatComplex operator-(const cuM2x2FloatComplex& a,
                                               const cuM2x2FloatComplex& b) {
  return cuM2x2FloatComplex(cuCsubf(a[0], b[0]), cuCsubf(a[1], b[1]),
                            cuCsubf(a[2], b[2]), cuCsubf(a[3], b[3]));
}

inline __device__ cuM2x2FloatComplex operator*(const cuM2x2FloatComplex& a,
                                               const cuM2x2FloatComplex& b) {
  return cuM2x2FloatComplex(cuCaddf(cuCmulf(a[0], b[0]), cuCmulf(a[1], b[2])),
                            cuCaddf(cuCmulf(a[0], b[1]), cuCmulf(a[1], b[3])),
                            cuCaddf(cuCmulf(a[2], b[0]), cuCmulf(a[3], b[2])),
                            cuCaddf(cuCmulf(a[2], b[1]), cuCmulf(a[3], b[3])));
}

inline __device__ cuM2x2DoubleComplex operator*(const cuM2x2DoubleComplex& a,
                                                const cuM2x2DoubleComplex& b) {
  return cuM2x2DoubleComplex(cuCadd(cuCmul(a[0], b[0]), cuCmul(a[1], b[2])),
                             cuCadd(cuCmul(a[0], b[1]), cuCmul(a[1], b[3])),
                             cuCadd(cuCmul(a[2], b[0]), cuCmul(a[3], b[2])),
                             cuCadd(cuCmul(a[2], b[1]), cuCmul(a[3], b[3])));
}

inline __device__ cuM2x2FloatComplex cuConj(const cuM2x2FloatComplex& x) {
  return cuM2x2FloatComplex(cuConjf(x[0]), cuConjf(x[2]), cuConjf(x[1]),
                            cuConjf(x[3]));
}

inline __device__ cuM2x2DoubleComplex cuConj(const cuM2x2DoubleComplex x) {
  return cuM2x2DoubleComplex(cuConj(x[0]), cuConj(x[2]), cuConj(x[1]),
                             cuConj(x[3]));
}

inline __device__ cuM2x2DoubleComplex
make_cuM2x2ComplexDouble(const cuM2x2FloatComplex& x) {
  return cuM2x2DoubleComplex(
      make_cuDoubleComplex(x[0]), make_cuDoubleComplex(x[1]),
      make_cuDoubleComplex(x[2]), make_cuDoubleComplex(x[3])

  );
}

using cuM2x2FloatComplexDiagonal = cuM2x2Diagonal<cuFloatComplex>;
using cuM2x2DoubleComplexDiagonal = cuM2x2Diagonal<cuDoubleComplex>;

inline __device__ cuM2x2FloatComplex
operator*(const cuM2x2FloatComplexDiagonal& a, const cuM2x2FloatComplex& b) {
  return cuM2x2FloatComplex(cuCmulf(a[0], b[0]), cuCmulf(a[0], b[1]),
                            cuCmulf(a[1], b[2]), cuCmulf(a[1], b[3]));
}

inline __device__ cuM2x2DoubleComplex
operator*(const cuM2x2DoubleComplexDiagonal& a, const cuM2x2DoubleComplex& b) {
  return cuM2x2DoubleComplex(cuCmul(a[0], b[0]), cuCmul(a[0], b[1]),
                             cuCmul(a[1], b[2]), cuCmul(a[1], b[3]));
}

#endif  // DP3_DDECAL_GAIN_SOLVERS_KERNELS_MATRIXCOMPLEX2X2_H_
