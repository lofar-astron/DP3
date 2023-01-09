// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
#define AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H

/**
 * @file Implements a Diagonal 2x2 Matrix with complex float values.
 *
 * This class is based on @ref aocommon::MC2x2FDiag but uses AVX-128
 * instructions.
 *
 * @note The class only implements a subset of the matrix operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon when the class has matured further.
 */

#include "common/avx256/VectorComplexFloat2.h"

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class DiagonalMatrixComplexFloat2x2 {
 public:
  [[nodiscard]] DiagonalMatrixComplexFloat2x2() noexcept = default;

  [[nodiscard]] /* implicit */ DiagonalMatrixComplexFloat2x2(
      VectorComplexFloat2 data) noexcept
      : data_{data} {}

  [[nodiscard]] explicit DiagonalMatrixComplexFloat2x2(
      const std::complex<float> a, const std::complex<float> b) noexcept
      : data_{a, b} {}

  [[nodiscard]] explicit DiagonalMatrixComplexFloat2x2(
      const std::complex<float> matrix[2]) noexcept
      : data_{VectorComplexFloat2{std::addressof(matrix[0])}} {}

  [[nodiscard]] std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] explicit operator __m128() const noexcept {
    return static_cast<__m128>(data_);
  }

  [[nodiscard]] DiagonalMatrixComplexFloat2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  [[nodiscard]] DiagonalMatrixComplexFloat2x2 HermitianTranspose()
      const noexcept {
    // The transpose has no effect for a diagonal matrix.
    return Conjugate();
  }

  DiagonalMatrixComplexFloat2x2 operator+=(
      DiagonalMatrixComplexFloat2x2 value) noexcept {
    return data_ += value.data_;
  }

  [[nodiscard]] friend bool operator==(
      DiagonalMatrixComplexFloat2x2 lhs,
      DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  DiagonalMatrixComplexFloat2x2 value) {
    output << "[{" << value[0] << ", " << std::complex<float>{} << "}, {"
           << std::complex<float>{} << ", " << value[1] << "}]";
    return output;
  }

  //
  // Deprecated operations
  //
  // The are resembling operations but use names not conforming to Google
  // Style or use named operations instead of operator overloading.
  //

  // RAP-133 enabled diagnostic [[deprecated("Use HermitianTranspose")]]
  [[nodiscard]] DiagonalMatrixComplexFloat2x2 HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  VectorComplexFloat2 data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
