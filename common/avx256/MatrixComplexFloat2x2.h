// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2X2_H
#define AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2X2_H

/**
 * @file Implements a 2x2 Matrix with complex float values.
 *
 * This class is based on @ref aocommon::MC2x2F but uses AVX-256 instructions.
 *
 * @note The class only implements a subset of the matrix operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon when the class has matured further.
 */

#include "common/avx256/DiagonalMatrixComplexFloat2x2.h"
#include "common/avx256/VectorComplexFloat4.h"

#include <aocommon/matrix2x2.h>

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class MatrixComplexFloat2x2 {
 public:
  [[nodiscard]] MatrixComplexFloat2x2() noexcept = default;

  [[nodiscard]] /* implicit */ MatrixComplexFloat2x2(
      VectorComplexFloat4 data) noexcept
      : data_{data} {}

  [[nodiscard]] explicit MatrixComplexFloat2x2(std::complex<float> a,
                                               std::complex<float> b,
                                               std::complex<float> c,
                                               std::complex<float> d) noexcept
      : data_{a, b, c, d} {}

  [[nodiscard]] explicit MatrixComplexFloat2x2(
      const std::complex<float> matrix[4]) noexcept
      : data_(matrix) {}

  [[nodiscard]] explicit MatrixComplexFloat2x2(
      const std::complex<double> matrix[4]) noexcept
      : data_(matrix) {}

  [[nodiscard]] explicit MatrixComplexFloat2x2(
      const aocommon::MC2x2F& matrix) noexcept
      : data_(matrix.Data()) {}

  [[nodiscard]] std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] explicit operator MC2x2F() const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return {data_[0], data_[1], data_[2], data_[3]};
  }

  [[nodiscard]] MatrixComplexFloat2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  [[nodiscard]] MatrixComplexFloat2x2 Transpose() const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return MatrixComplexFloat2x2{data_[0], data_[2], data_[1], data_[3]};
  }

  [[nodiscard]] MatrixComplexFloat2x2 HermitianTranspose() const noexcept {
    return Transpose().Conjugate();
  }

  explicit operator __m256() const noexcept {
    return static_cast<__m256>(data_);
  }

  /// @returns the Frobenius norm of the matrix.
  [[nodiscard]] float Norm() const noexcept {
    // This uses the same basic idea as MatrixComplexDouble2x2::Norm except
    // that the underlaying data is stored in one __m256d value.

    // Note this function seems slower than expected.
    // MatrixComplexDouble2x2::Norm is faster than this function. It is
    // still faster than the scalare version. It would be nice to improve
    // this in the future.
    //
    //
    __m256 tmp = static_cast<__m256>(data_);
    tmp *= tmp;

    // For the addition we deviate from the double version.
    // | a     | b     | c     | d     |
    // | e     | f     | g     | h     |
    // ---------------------------------+
    // | a + e | b + f | c + g | d + h |

    __m128 ret = _mm256_castps256_ps128(tmp);
    ret += _mm256_extractf128_ps(tmp, 1);

    // | a + e         | b + f         |
    // | c + g         | d + h         |
    // ---------------------------------+
    // | a + e + c + g | b + f + d + h |

    // ret = _mm_hadd_ps(ret, ret);

    ret += _mm_movehl_ps(ret, ret);

    // | a + e + c + g                 |
    // | b + f + d + h                 |
    // ---------------------------------+
    // | a + e + c + g + b + f + d + h |

    return ret[0] + ret[1];
  }

  /** Returns the sum of the diagonal elements. */
  [[nodiscard]] std::complex<float> Trace() const noexcept {
    // Trace = M[0] + M[3]

    __m256 ret = static_cast<__m256>(data_);
    ret +=
        // Moves M[3] to the location of M[0] and adds it to ret.
        _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(ret), 0b11));
    return {ret[0], ret[1]};
  }

  /** Assign data stored by 2x2 complex matrix to destination buffer */
  void AssignTo(std::complex<float>* destination) const noexcept {
    data_.AssignTo(destination);
  }

  [[nodiscard]] static MatrixComplexFloat2x2 Unity() noexcept {
    return MatrixComplexFloat2x2{
        std::complex<float>(1.0f, 0.0f), std::complex<float>(0.0f, 0.0f),
        std::complex<float>(0.0f, 0.0f), std::complex<float>(1.0f, 0.0f)};
  }

  [[nodiscard]] static MatrixComplexFloat2x2 NaN() noexcept {
    return MatrixComplexFloat2x2{
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()}};
  }

  MatrixComplexFloat2x2& operator+=(MatrixComplexFloat2x2 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  MatrixComplexFloat2x2& operator-=(MatrixComplexFloat2x2 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator+(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator-(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return lhs -= rhs;
  }

  MatrixComplexFloat2x2& operator*=(MatrixComplexFloat2x2 value) noexcept {
    data_ = data_ * value.data_;
    return *this;
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
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

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, std::complex<float> rhs) noexcept {
    return lhs.data_ * rhs;
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      std::complex<float> lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] DiagonalMatrixComplexFloat2x2 Diagonal() const noexcept {
    return VectorComplexFloat2{data_[0], data_[3]};
  }

  [[nodiscard]] bool Invert() noexcept {
    // Note it would be possible to optimize further.
    VectorComplexFloat2 a{data_[0], data_[1]};
    VectorComplexFloat2 b{data_[3], data_[2]};
    VectorComplexFloat2 c = a * b;

    std::complex<float> d = c[0] - c[1];
    if (d == std::complex<float>{}) return false;

    float n = std::norm(d);
    d.imag(-d.imag());
    __m256 reprocical = _mm256_setr_ps(d.real(), d.imag(), d.real(), d.imag(),
                                       d.real(), d.imag(), d.real(), d.imag());
    reprocical = _mm256_div_ps(reprocical, _mm256_set1_ps(n));

    // std::swap(data[0],data[3]);
    // Using the fact that extracting as a double, the value has the number of
    // bits is the same as a complex gloat
    __m256d data = _mm256_castps_pd(static_cast<__m256>(data_));
    data = _mm256_permute4x64_pd(data, 0b00'10'01'11);
    __m256 result = _mm256_castpd_ps(data);

    // data[0] = data[0]
    // data[1] = -data[1]
    // data[2] = -data[2]
    // data[3] = data[3]
    __m256 mask = _mm256_setr_ps(0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0);
    result = _mm256_xor_ps(result, mask);

    data_ = VectorComplexFloat4{result} * VectorComplexFloat4{reprocical};

    return true;
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    // basically this does the same as
    // return lhs.data_ * VectorComplexFloat4{rhs[0], rhs[0], rhs[1], rhs[1]};
    // but this gives better codegen and is slightly faster.

    // Since there is no direct way to interleave 2 32-bit float pairs, use the
    // method to interleave 1 64-bit value. (The values are unchanged.)

    // __m128 lo {rhs[0], rhs[0]}
    __m128 lo =
        _mm_castpd_ps(_mm_unpacklo_pd(_mm_castps_pd(static_cast<__m128>(rhs)),
                                      _mm_castps_pd(static_cast<__m128>(rhs))));
    // __m128 hi {rhs[1], rhs[1]}
    __m128 hi =
        _mm_castpd_ps(_mm_unpackhi_pd(_mm_castps_pd(static_cast<__m128>(rhs)),
                                      _mm_castps_pd(static_cast<__m128>(rhs))));
    return lhs.data_ * VectorComplexFloat4{_mm256_set_m128(hi, lo)};
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      DiagonalMatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] friend bool operator==(MatrixComplexFloat2x2 lhs,
                                       MatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  MatrixComplexFloat2x2 value) {
    output << "[{" << value[0] << ", " << value[1] << "}, {" << value[2] << ", "
           << value[3] << "}]";
    return output;
  }

 private:
  VectorComplexFloat4 data_;
};

/// MC2x2Base compatibility wrapper.
inline MatrixComplexFloat2x2 HermTranspose(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.HermitianTranspose();
}

/**
 * MC2x2Base compatibility wrapper.
 *
 * @returns the Frobenius norm of the matrix.
 */
inline float Norm(MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Norm();
}

/** Returns the sum of the diagonal elements. */
inline std::complex<float> Trace(MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Trace();
}

[[nodiscard]] inline DiagonalMatrixComplexFloat2x2 Diagonal(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Diagonal();
}

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2X2_H
