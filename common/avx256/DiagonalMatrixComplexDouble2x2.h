// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
#define AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H

/**
 * @file Implements a Diagonal 2x2 Matrix with complex double values.
 *
 * This class is based on @ref aocommon::MC2x2Diag but uses AVX-256
 * instructions.
 *
 * @warning All functions in this header need to use a target attribute
 * like @c [[gnu::target("avx2,fma")]]. When this is not done the GCC
 * doesn't adhere to the proper ABI leading to broken code.
 *
 * @note The class only implements a subset of the matrix operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon when the class has matured further.
 */

#include "common/avx256/VectorComplexDouble2.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>

namespace aocommon::Avx256 {

class DiagonalMatrixComplexDouble2x2 {
 public:
  [[nodiscard]] [[gnu::target(
      "avx2,fma")]] DiagonalMatrixComplexDouble2x2() noexcept = default;

  [[nodiscard]] [[gnu::target("avx2,fma")]] /* implicit */
  DiagonalMatrixComplexDouble2x2(VectorComplexDouble2 data) noexcept
      : data_{data} {}

  [[nodiscard]] [[gnu::target(
      "avx2,fma")]] explicit DiagonalMatrixComplexDouble2x2(const std::
                                                                complex<double>
                                                                    a,
                                                            const std::complex<
                                                                double>
                                                                b) noexcept
      : data_{a, b} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] VectorComplexDouble2 Data()
      const noexcept {
    return data_;
  }

  [[nodiscard]] [[gnu::target(
      "avx2,fma")]] explicit DiagonalMatrixComplexDouble2x2(const std::
                                                                complex<double>
                                                                    matrix
                                                                        [2]) noexcept
      : data_{VectorComplexDouble2{std::addressof(matrix[0])}} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] std::complex<double> operator[](
      size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] DiagonalMatrixComplexDouble2x2
  Conjugate() const noexcept {
    return data_.Conjugate();
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] DiagonalMatrixComplexDouble2x2
  HermitianTranspose() const noexcept {
    // The transpose has no effect for a diagonal matrix.
    return Conjugate();
  }

  [[gnu::target("avx2,fma")]] DiagonalMatrixComplexDouble2x2 operator+=(
      DiagonalMatrixComplexDouble2x2 value) noexcept {
    return data_ += value.data_;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      DiagonalMatrixComplexDouble2x2 lhs,
      DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  [[gnu::target("avx2,fma")]] friend std::ostream& operator<<(
      std::ostream& output, DiagonalMatrixComplexDouble2x2 value) {
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
  [[nodiscard]] [[gnu::target("avx2,fma")]] DiagonalMatrixComplexDouble2x2
  HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  VectorComplexDouble2 data_;
};

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
