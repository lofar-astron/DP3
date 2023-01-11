// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
#define AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H

/**
 * @file Implements a Diagnoal 2x2 Matrix with complex double values.
 *
 * This class is based on @ref aocommon::MC2x2Diag but uses AVX-256
 * instructions.
 *
 * @note The class only implements a subset of the matrix operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */

#include "common/avx256/VectorComplexDouble2.h"

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class DiagonalMatrixComplexDouble2x2 {
 public:
  [[nodiscard]] DiagonalMatrixComplexDouble2x2() noexcept = default;

  [[nodiscard]] /* implicit */ DiagonalMatrixComplexDouble2x2(
      VectorComplexDouble2 data) noexcept
      : data_{data} {}

  [[nodiscard]] explicit DiagonalMatrixComplexDouble2x2(
      const std::complex<double> a, const std::complex<double> b) noexcept
      : data_{a, b} {}

  [[nodiscard]] explicit DiagonalMatrixComplexDouble2x2(
      const std::complex<double> matrix[2]) noexcept
      : data_{VectorComplexDouble2{std::addressof(matrix[0])}} {}

  [[nodiscard]] std::complex<double> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] DiagonalMatrixComplexDouble2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  [[nodiscard]] DiagonalMatrixComplexDouble2x2 HermitianTranspose()
      const noexcept {
    // The transpose has no effect for a diagonal matrix.
    return Conjugate();
  }

  [[nodiscard]] friend bool operator==(
      DiagonalMatrixComplexDouble2x2 lhs,
      DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  DiagonalMatrixComplexDouble2x2 value) {
    output << "[{" << value[0] << ", " << std::complex<double>{} << "}, {"
           << std::complex<double>{} << ", " << value[1] << "}]";
    return output;
  }

  //
  // Deprecated operations
  //
  // The are resembling operations but use names not conforming to Google
  // Style or use named operations instead of operator overloading.
  //

  // RAP-133 enabled diagnostic [[deprecated("Use HermitianTranspose")]]
  [[nodiscard]] DiagonalMatrixComplexDouble2x2 HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  VectorComplexDouble2 data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
