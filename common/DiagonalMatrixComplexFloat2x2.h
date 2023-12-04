// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
#define AOCOMMON_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H

#include "common/scalar/DiagonalMatrixComplexFloat2x2.h"
#include "common/avx256/DiagonalMatrixComplexFloat2x2.h"
#include "common/MatrixComplexFloat2x2.h"

#include <functional>

/**
 * @file AVX or Scalar dispatch code for the complex float diagonal matrix.
 *
 * The dispatching used in this file is similar to the dispatching in
 * @ref MatrixComplexFloat2x2.h. That file contains the documenation for the
 * design.
 */

namespace aocommon {

#ifndef __AVX2__
// AVX isn't present at build time.
using DiagonalMatrixComplexFloat2x2 = Scalar::DiagonalMatrixComplexFloat2x2;

#elif !defined(PORTABLE_BUILD)  // __AVX2__

// AVX is present at build time and non-portable builds.
using DiagonalMatrixComplexFloat2x2 = Avx256::DiagonalMatrixComplexFloat2x2;

#else  // __AVX2__

// AVX is present at build time and non-portable builds.
class DiagonalMatrixComplexFloat2x2 {
  enum class Dispatch { Scalar, Avx256 };

  [[nodiscard]] [[gnu::target("default")]] static Dispatch GetDispatch() {
    return Dispatch::Scalar;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] static Dispatch GetDispatch() {
    return Dispatch::Avx256;
  }

  /// Converting constructor using a scalar matrix.
  DiagonalMatrixComplexFloat2x2(Scalar::DiagonalMatrixComplexFloat2x2 scalar)
      : scalar_{scalar} {}

  /// Converting constructor using an AVX  matrix.
  DiagonalMatrixComplexFloat2x2(Avx256::DiagonalMatrixComplexFloat2x2 avx)
      : avx_{avx} {}

  template <class Function, class... Args>
  [[nodiscard]] auto ExecuteCall(Function function, Args&&... args) const {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return function(scalar_, std::forward<Args>(args)...);
      case Dispatch::Avx256:
        return function(avx_, std::forward<Args>(args)...);
    }
    __builtin_unreachable();
  }

 public:
  [[gnu::target("default")]] explicit DiagonalMatrixComplexFloat2x2(
      const std::complex<double> matrix[2]) noexcept
      : scalar_{matrix} {}

 private:
  union {
    Scalar::DiagonalMatrixComplexFloat2x2 scalar_;
    Avx256::DiagonalMatrixComplexFloat2x2 avx_;
  };
};

#endif  // __AVX2__

}  // namespace aocommon

#endif  // AOCOMMON_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
