// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMPLEX_H_
#define DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMPLEX_H_

#include <cuComplex.h>

__host__ __device__ static __inline__ cuDoubleComplex make_cuDoubleComplex(
    const cuFloatComplex& a) {
  return make_cuDoubleComplex(a.x, a.y);
}

/*
 * Taken the below utility functions from
 * https://forums.developer.nvidia.com/t/additional-cucomplex-functions-cucnorm-cucsqrt-cucexp-and-some-complex-double-functions/36892
 */
__host__ __device__ static __inline__ cuDoubleComplex cuCadd(cuDoubleComplex x,
                                                             double y) {
  return make_cuDoubleComplex(cuCreal(x) + y, cuCimag(x));
}
__host__ __device__ static __inline__ cuDoubleComplex cuCdiv(cuDoubleComplex x,
                                                             double y) {
  return make_cuDoubleComplex(cuCreal(x) / y, cuCimag(x) / y);
}
__host__ __device__ static __inline__ cuDoubleComplex cuCmul(cuDoubleComplex x,
                                                             double y) {
  return make_cuDoubleComplex(cuCreal(x) * y, cuCimag(x) * y);
}
__host__ __device__ static __inline__ cuDoubleComplex cuCsub(cuDoubleComplex x,
                                                             double y) {
  return make_cuDoubleComplex(cuCreal(x) - y, cuCimag(x));
}

__host__ __device__ static __inline__ cuDoubleComplex cuCexp(
    cuDoubleComplex x) {
  double factor = exp(x.x);
  return make_cuDoubleComplex(factor * cos(x.y), factor * sin(x.y));
}

/*
 * Cuda complex implementation of std::arg
 * https://en.cppreference.com/w/cpp/numeric/complex/arg
 */
__device__ static __inline__ float cuCarg(const cuDoubleComplex& z) {
  return atan2(cuCimag(z), cuCreal(z));
}

/*
 * Cuda complex implementation of std::polar
 * https://en.cppreference.com/w/cpp/numeric/complex/polar
 */
__device__ static __inline__ cuDoubleComplex cuCpolar(const double r,
                                                      const double z) {
  return make_cuDoubleComplex(r * cos(z), r * sin(z));
}

template <typename T>
__device__ static __inline__ double cuNorm(const T& a) {
  return a.x * a.x + a.y * a.y;
}

#endif  // DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMPLEX_H_