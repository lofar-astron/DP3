// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_MATRIX_COMPLEX_DOUBLE_2X2_H
#define AOCOMMON_AVX256_MATRIX_COMPLEX_DOUBLE_2X2_H

/**
 * @file Implements a 2x2 Matrix with complex double values.
 *
 * This class is based on @ref aocommon::MC2x2 but uses AVX-256 instructions.
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

#include "common/avx256/DiagonalMatrixComplexDouble2x2.h"
#include "common/avx256/MatrixComplexFloat2x2.h"
#include "common/avx256/VectorComplexDouble2.h"

#include <aocommon/matrix2x2.h>

#include <array>
#include <cassert>
#include <complex>
#include <immintrin.h>
#include <limits>
#include <ostream>

namespace aocommon::Avx256 {

class MatrixComplexDouble2x2 {
 public:
  [[nodiscard]] [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2() noexcept =
      default;

  [[nodiscard]] [[gnu::target("avx2,fma")]] /* implicit */
  MatrixComplexDouble2x2(std::array<VectorComplexDouble2, 2> data) noexcept
      : data_{data} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit MatrixComplexDouble2x2(
      std::complex<double> a, std::complex<double> b, std::complex<double> c,
      std::complex<double> d) noexcept
      : data_{{VectorComplexDouble2{a, b}, VectorComplexDouble2{c, d}}} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit MatrixComplexDouble2x2(
      const std::complex<float> matrix[4]) noexcept
      : data_{{VectorComplexDouble2{std::addressof(matrix[0])},
               VectorComplexDouble2{std::addressof(matrix[2])}}} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit MatrixComplexDouble2x2(
      MatrixComplexFloat2x2 matrix) noexcept {
    __m256 tmp = static_cast<__m256>(matrix);

    __m128 lo = _mm256_castps256_ps128(tmp);
    __m128 hi = _mm256_extractf128_ps(tmp, 1);
    data_[0] = VectorDouble4{_mm256_cvtps_pd(lo)};
    data_[1] = VectorDouble4{_mm256_cvtps_pd(hi)};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit MatrixComplexDouble2x2(
      const std::complex<double> matrix[4]) noexcept
      : data_{{VectorComplexDouble2{std::addressof(matrix[0])},
               VectorComplexDouble2{std::addressof(matrix[2])}}} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit MatrixComplexDouble2x2(
      const aocommon::MC2x2& matrix) noexcept
      : MatrixComplexDouble2x2(matrix.Data()) {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] std::complex<double> operator[](
      size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    size_t array = index / 2;
    index %= 2;
    return data_[array][index];
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2 Conjugate()
      const noexcept {
    return std::array<VectorComplexDouble2, 2>{data_[0].Conjugate(),
                                               data_[1].Conjugate()};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2 Transpose()
      const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return MatrixComplexDouble2x2{(*this)[0], (*this)[2], (*this)[1],
                                  (*this)[3]};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2
  HermitianTranspose() const noexcept {
    return Transpose().Conjugate();
  }

  /// @returns the Frobenius norm of the matrix.
  [[nodiscard]] [[gnu::target("avx2,fma")]] double Norm() const noexcept {
    // Norm Matrix Complex 2x2
    // Norm(a) + Norm(b) + Norm(c) + Norm(d)
    //
    // Norm is using C++'s definition of std::complex<T>::norm(). This norm is
    // also known as the 'field norm' or 'absolute square'. Norm is defined as
    // a.re * a.re + a.im * a.im
    //
    // Note if we want to do this accoring to the rules above some shuffing
    // needs to be done. Instead we can consider the underlaying data an array
    // of 8 doubles. Then Norm becomes
    //
    // -- 7
    // \.
    //  .  a[n] * a[n]
    // /
    // -- n = 0
    //
    // and no shuffling in needed instead use the following algorithm
    //
    // hi = data_[0]
    // lo = data_[1]
    //
    // hi = hi * hi
    // lo = lo * lo
    //
    // tmp = hi + lo
    // ret = std::accumulate(&tmp[0], &tmp[4], 0.0); // not possible in C++
    //
    // instead of calculating tmp as described it can be done by
    // hi = lo * lo + hi

    __m256d hi = static_cast<__m256d>(data_[0]);
    __m256d lo = static_cast<__m256d>(data_[1]);

    hi *= hi;
    hi = _mm256_fmadd_pd(lo, lo, hi);

    // Summing the 4 elements in hi can be simply done by
    // return hi[0] + hi[1] + hi[2] + hi[3]
    //
    // however this is slow, it's more efficient to permutate the data and use
    // vector adding. The instruction set has a hadd operation, but this is
    // slow too. Instead use register permutations and additons. The entries
    // marked with - in the table mean we don't care about the contents. The
    // result will be stored in hi[0]:
    //
    // hi | a             | b     | c | d |
    // lo | c             | d     | - | - |
    //    --------------------------------- +
    // hi | a + c         | b + d | - | - |
    // lo | b + d         | -     | - | - |
    //    --------------------------------- +
    // hi | a + c + b + d | -     | - | - |

    lo = _mm256_permute4x64_pd(hi, 0b11'10);
    hi += lo;

    __m128d ret = _mm256_castpd256_pd128(hi);
    ret += _mm_permute_pd(ret, 0b01);
    return ret[0];
  }

  /** Returns the sum of the diagonal elements. */
  [[nodiscard]] [[gnu::target("avx2,fma")]] std::complex<double> Trace()
      const noexcept {
    // Trace = M[0] + M[3]

    // Extract M[0] and M[1] as 128-bit AVX vector.
    __m128d ret = _mm256_castpd256_pd128(static_cast<__m256d>(data_[0]));
    // Extracts M[3] and M[4] as 128-bit AVX vector and adds it to ret.
    ret += _mm256_extractf128_pd(static_cast<__m256d>(data_[1]), 1);
    return {ret[0], ret[1]};
  }

  /** Assign data stored by 2x2 complex matrix to destination buffer */
  [[gnu::target("avx2,fma")]] void AssignTo(
      std::complex<double>* destination) const noexcept {
    data_[0].AssignTo(destination);
    destination += 2;
    data_[1].AssignTo(destination);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] static MatrixComplexDouble2x2
  Unity() noexcept {
    return MatrixComplexDouble2x2{
        std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0),
        std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] static MatrixComplexDouble2x2
  NaN() noexcept {
    return MatrixComplexDouble2x2{
        std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()},
        std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()},
        std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()},
        std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()}};
  }

  [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2& operator+=(
      MatrixComplexDouble2x2 value) noexcept {
    data_[0] += value.data_[0];
    data_[1] += value.data_[1];
    return *this;
  }

  [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2& operator-=(
      MatrixComplexDouble2x2 value) noexcept {
    data_[0] -= value.data_[0];
    data_[1] -= value.data_[1];
    return *this;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator+(MatrixComplexDouble2x2 lhs, MatrixComplexDouble2x2 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator-(MatrixComplexDouble2x2 lhs, MatrixComplexDouble2x2 rhs) noexcept {
    return lhs -= rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator*(MatrixComplexDouble2x2 lhs, MatrixComplexDouble2x2 rhs) noexcept {
    // The 2x2 matrix multiplication is done using the following algorithm.

    // High:
    // ret.a = lhs.a * rhs.a + lhs.b * rhs.c
    // ret.b = lhs.a * rhs.b + lhs.b * rhs.d
    //       | hc1   | hc2   | hc3   | hc4   |
    //       | hs1           | hs2           |

    // Low:
    // ret.c = lhs.c * rhs.a + lhs.d * rhs.c
    // ret.d = lhs.c * rhs.b + lhs.d * rhs.d
    //       | lc1   | lc2   | lc3   | lc4   |
    //       | ls1           | ls2           |

    // High:
    VectorComplexDouble2 hc1{lhs[0], lhs[0]};
    VectorComplexDouble2 hc2{rhs[0], rhs[1]};
    VectorComplexDouble2 hs1 = hc1 * hc2;

    VectorComplexDouble2 hc3{lhs[1], lhs[1]};
    VectorComplexDouble2 hc4{rhs[2], rhs[3]};
    VectorComplexDouble2 hs2 = hc3 * hc4;

    VectorComplexDouble2 hr = hs1 + hs2;

    // Low:
    VectorComplexDouble2 lc1{lhs[2], lhs[2]};
    VectorComplexDouble2 lc2{rhs[0], rhs[1]};
    VectorComplexDouble2 ls1 = lc1 * lc2;

    VectorComplexDouble2 lc3{lhs[3], lhs[3]};
    VectorComplexDouble2 lc4{rhs[2], rhs[3]};
    VectorComplexDouble2 ls2 = lc3 * lc4;

    VectorComplexDouble2 lr = ls1 + ls2;

    return std::array<VectorComplexDouble2, 2>{hr, lr};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] DiagonalMatrixComplexDouble2x2
  Diagonal() const noexcept {
    return DiagonalMatrixComplexDouble2x2(data_[0][0], data_[1][1]);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator*(MatrixComplexDouble2x2 lhs, std::complex<double> rhs) noexcept {
    return std::array<VectorComplexDouble2, 2>{lhs.data_[0] * rhs,
                                               lhs.data_[1] * rhs};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator*(std::complex<double> lhs, MatrixComplexDouble2x2 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator*(MatrixComplexDouble2x2 lhs,
            DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return std::array<VectorComplexDouble2, 2>{lhs.data_[0] * rhs[0],
                                               lhs.data_[1] * rhs[1]};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend MatrixComplexDouble2x2
  operator*(DiagonalMatrixComplexDouble2x2 lhs,
            MatrixComplexDouble2x2 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      MatrixComplexDouble2x2 lhs, MatrixComplexDouble2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2& operator=(
      const aocommon::MC2x2& matrix) noexcept {
    *this = MatrixComplexDouble2x2(matrix.Data());
    return *this;
  }

  [[gnu::target("avx2,fma")]] friend std::ostream& operator<<(
      std::ostream& output, MatrixComplexDouble2x2 value) {
    output << "[{" << value[0] << ", " << value[1] << "}, {" << value[2] << ", "
           << value[3] << "}]";
    return output;
  }

  // RAP-133 enabled diagnostic [[deprecated("Use HermitianTranspose")]]
  [[nodiscard]] [[gnu::target("avx2,fma")]] MatrixComplexDouble2x2
  HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  std::array<VectorComplexDouble2, 2> data_;
};

// RAP-133 enabled diagnostic
// [[deprecated("Use MatrixComplexDouble2x2::HermitianTranspose")]]
/// MC2x2Base compatibility wrapper.
[[nodiscard]] [[gnu::target("avx2,fma")]] inline MatrixComplexDouble2x2
HermTranspose(MatrixComplexDouble2x2 matrix) noexcept {
  return matrix.HermitianTranspose();
}

/**
 * MC2x2Base compatibility wrapper.
 *
 * @returns the Frobenius norm of the matrix.
 */
[[nodiscard]] [[gnu::target("avx2,fma")]] inline double Norm(
    MatrixComplexDouble2x2 matrix) noexcept {
  return matrix.Norm();
}

/** Returns the sum of the diagonal elements. */
[[nodiscard]] [[gnu::target("avx2,fma")]] inline std::complex<double> Trace(
    MatrixComplexDouble2x2 matrix) noexcept {
  return matrix.Trace();
}

[[nodiscard]] [[gnu::target("avx2,fma")]] inline DiagonalMatrixComplexDouble2x2
Diagonal(MatrixComplexDouble2x2 matrix) noexcept {
  return matrix.Diagonal();
}

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_MATRIX_COMPLEX_DOUBLE_2X2_H
