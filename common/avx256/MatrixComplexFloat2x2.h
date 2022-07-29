// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2_X_2_H
#define AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2_X_2_H

/**
 * @file Implements a 2x2 Matrix with complex float values.
 *
 * This class is based on @ref aocommon::MC2x2F but uses AVX-256 instructions.
 *
 * @note The class only implements a subset of the matrix operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */

#include "common/avx256/VectorComplexFloat4.h"

#include <aocommon/matrix2x2.h>

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class MaxtrixComplexFloat2x2 {
 public:
  /* implicit */ MaxtrixComplexFloat2x2(VectorComplexFloat4 data) noexcept
      : data_{data} {}

  explicit MaxtrixComplexFloat2x2(std::complex<float> a, std::complex<float> b,
                                  std::complex<float> c,
                                  std::complex<float> d) noexcept
      : data_{a, b, c, d} {}

  explicit MaxtrixComplexFloat2x2(const std::complex<float> matrix[4]) noexcept
      : data_(matrix) {}

  explicit MaxtrixComplexFloat2x2(const std::complex<double> matrix[4]) noexcept
      : data_(matrix) {}

  explicit MaxtrixComplexFloat2x2(const aocommon::MC2x2F& matrix) noexcept
      : data_(matrix.Data()) {}

  [[nodiscard]] std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] explicit operator MC2x2F() const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return {data_[0], data_[1], data_[2], data_[3]};
  }

  [[nodiscard]] MaxtrixComplexFloat2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  [[nodiscard]] MaxtrixComplexFloat2x2 Transpose() const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return MaxtrixComplexFloat2x2{data_[0], data_[2], data_[1], data_[3]};
  }

  [[nodiscard]] MaxtrixComplexFloat2x2 HermitianTranspose() const noexcept {
    return Transpose().Conjugate();
  }

  [[nodiscard]] friend MaxtrixComplexFloat2x2 operator*(
      MaxtrixComplexFloat2x2 lhs, MaxtrixComplexFloat2x2 rhs) noexcept {
    // The 2x2 matrix multiplication is done using the following algorithm.
    // ret.a = lhs.a * rhs.a + lhs.b * rhs.c
    // ret.b = lhs.a * rhs.b + lhs.b * rhs.d
    // ret.c = lhs.c * rhs.a + lhs.d * rhs.c
    // ret.d = lhs.c * rhs.b + lhs.d * rhs.d
    //       | c1    | c2    | c3    | c4    |
    //       | s1            | s2            |
    //

    VectorComplexFloat4 c1{lhs[0], lhs[0], lhs[2], lhs[2]};
    VectorComplexFloat4 c2{rhs[0], rhs[1], rhs[0], rhs[1]};
    VectorComplexFloat4 s1 = c1 * c2;

    VectorComplexFloat4 c3{lhs[1], lhs[1], lhs[3], lhs[3]};
    VectorComplexFloat4 c4{rhs[2], rhs[3], rhs[2], rhs[3]};
    VectorComplexFloat4 s2 = c3 * c4;

    return s1 + s2;
  }

  [[nodiscard]] friend bool operator==(MaxtrixComplexFloat2x2 lhs,
                                       MaxtrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  MaxtrixComplexFloat2x2 value) {
    output << "[{" << value[0] << ", " << value[1] << "}, {" << value[2] << ", "
           << value[3] << "}]";
    return output;
  }

 private:
  VectorComplexFloat4 data_;
};

/// MC2x2Base compatibility wrapper.
inline MaxtrixComplexFloat2x2 HermTranspose(
    MaxtrixComplexFloat2x2 matrix) noexcept {
  return matrix.HermitianTranspose();
}

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2_X_2_H