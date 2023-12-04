// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
#define AOCOMMON_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H

#include "common/scalar/DiagonalMatrixComplexDouble2x2.h"
#ifdef __AVX2__
#include "common/avx256/DiagonalMatrixComplexDouble2x2.h"
#endif
#include "common/MatrixComplexDouble2x2.h"

#include <functional>

/**
 * @file AVX or Scalar dispatch code for the complex double diagonal matrix.
 *
 * The dispatching used in this file is similar to the dispatching in
 * @ref MatrixComplexDouble2x2.h. That file contains the documenation for the
 * design.
 */

namespace aocommon {

#ifndef __AVX2__
// AVX isn't present at build time.
using DiagonalMatrixComplexDouble2x2 = Scalar::DiagonalMatrixComplexDouble2x2;

#elif !defined(PORTABLE_BUILD)  // __AVX2__

// AVX is present at build time and non-portable builds.
using DiagonalMatrixComplexDouble2x2 = Avx256::DiagonalMatrixComplexDouble2x2;

#else  // __AVX2__

// AVX is present at build time and non-portable builds.
class DiagonalMatrixComplexDouble2x2 {
  enum class Dispatch { Scalar, Avx256 };

  [[nodiscard]] [[gnu::target("default")]] static Dispatch GetDispatch() {
    return Dispatch::Scalar;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] static Dispatch GetDispatch() {
    return Dispatch::Avx256;
  }

  /// Converting constructor using a scalar matrix.
  DiagonalMatrixComplexDouble2x2(Scalar::DiagonalMatrixComplexDouble2x2 scalar)
      : scalar_{scalar} {}

  /// Converting constructor using an AVX  matrix.
  DiagonalMatrixComplexDouble2x2(Avx256::DiagonalMatrixComplexDouble2x2 avx)
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
  [[gnu::target("default")]] explicit DiagonalMatrixComplexDouble2x2(
      const std::complex<double> matrix[2]) noexcept
      : scalar_{matrix} {}

  //
  // Operations
  //
  // The use the private helper functions to dispatch the call to the proper
  // underlying type.
  //

  [[nodiscard]] std::complex<double> operator[](size_t index) const noexcept {
    return ExecuteCall(
        [](const auto& object, size_t index) { return object[index]; }, index);
  }

  [[nodiscard]] DiagonalMatrixComplexDouble2x2 HermitianTranspose()
      const noexcept {
    return ExecuteCall([](auto object) -> DiagonalMatrixComplexDouble2x2 {
      if constexpr (std::is_same_v<decltype(object),
                                   Scalar::DiagonalMatrixComplexDouble2x2>)
        return object.HermTranspose();
      else
        return object.HermitianTranspose();
    });
  }

  //
  // Deprecated operations
  //
  // The are resembling operations but use names not conforming to Google
  // Style or use named operations instead of operator overloading.
  //

  [[deprecated(
      "Use HermitianTranspose")]] [[nodiscard]] DiagonalMatrixComplexDouble2x2
  HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  union {
    Scalar::DiagonalMatrixComplexDouble2x2 scalar_;
    Avx256::DiagonalMatrixComplexDouble2x2 avx_;
  };
};

#endif  // __AVX2__

}  // namespace aocommon

#endif  // AOCOMMON_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
