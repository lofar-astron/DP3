// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H
#define AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H

/**
 * @file Implements a Vector with 4 complex float values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::MatrixComplexFloat2x2, but can be used by itself.
 *
 * @warning All functions in this header need to use a target attribute
 * like @c [[gnu::target("avx2,fma")]]. When this is not done the GCC
 * doesn't adhere to the proper ABI leading to broken code.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon when the class has matured further.
 */

#include "common/avx256/VectorFloat8.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>

namespace aocommon::Avx256 {

class VectorComplexFloat4 {
 public:
  [[gnu::target("avx2,fma")]] VectorComplexFloat4() noexcept
      : data_{_mm256_setzero_ps()} {}

  [[gnu::target("avx2,fma")]] VectorComplexFloat4(VectorFloat8 data) noexcept
      : data_{data} {}

  [[gnu::target("avx2,fma")]] explicit VectorComplexFloat4(
      std::complex<float> a, std::complex<float> b, std::complex<float> c,
      std::complex<float> d) noexcept
      : data_{VectorFloat8{a.real(), a.imag(), b.real(), b.imag(), c.real(),
                           c.imag(), d.real(), d.imag()}} {}

  [[gnu::target("avx2,fma")]] explicit VectorComplexFloat4(
      const std::complex<float> vector[4]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat8{
            reinterpret_cast<const float*>(std::addressof(vector[0]))}} {}

  [[gnu::target("avx2,fma")]] explicit VectorComplexFloat4(
      const std::complex<double> vector[4]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat8{
            reinterpret_cast<const double*>(std::addressof(vector[0]))}} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] std::complex<float> operator[](
      size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return {data_[2 * index], data_[2 * index + 1]};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit operator __m256()
      const noexcept {
    return data_.Value();
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] VectorComplexFloat4 Conjugate()
      const noexcept {
    // Xor-ing a float with  0.0 will not change the value.
    // Xor-ing a float with -0.0 will change the sign of the value.
    __m256 mask = _mm256_setr_ps(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0);
    return data_ ^ mask;
  }

  /** Assign data stored by 4 element complex vector to destination buffer. */
  [[gnu::target("avx2,fma")]] void AssignTo(
      std::complex<float>* destination) const noexcept {
    data_.AssignTo(reinterpret_cast<float*>(destination));
  }

  [[gnu::target("avx2,fma")]] void AssignTo(
      std::complex<double>* destination) const noexcept {
    data_.AssignTo(reinterpret_cast<double*>(destination));
  }

  [[gnu::target("avx2,fma")]] VectorComplexFloat4& operator+=(
      VectorComplexFloat4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  [[gnu::target("avx2,fma")]] VectorComplexFloat4& operator-=(
      VectorComplexFloat4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorComplexFloat4
  operator+(VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorComplexFloat4
  operator-(VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
    return lhs -= rhs;
  }

  /// Multiplies the elements of 2 vectors on parallel.
  ///
  /// r[0] = lhs[0] * rhs[0]
  /// r[1] = lhs[1] * rhs[1]
  /// r[2] = lhs[2] * rhs[2]
  /// r[3] = lhs[3] * rhs[3]
  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorComplexFloat4
  operator*(VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
    // The complex multiplication L * R is performed by:
    // L * R = (L Re * R Re) - (L Im * R Im) +
    //        ((L Re * R Im) + (L Im * R Re)) i

    // Transformed for AVX this becomes
    //
    // mul 1:
    //
    // Lv1 | L Re        | L Re        |
    // Lv2 | R Re        | R Im        | Note this is rhs
    //     ------------- * ----------- *
    // Lv3 | L Re * R Re | L Re * R Im |
    //
    // mul 2:
    //
    // Rv1 | L Im        | L Im        |
    // Rv2 | R Im        | R Re        |
    //     ------------- * ----------- *
    // Rv3 | L Im * R Im | L Im * R Re |
    //
    // add sub
    // Lv3 | L Re * R Re               | L Re * R Im               |
    // Rv3 | L Im * R Im               | L Im * R Re               |
    //     --------------------------- - ------------------------- +
    //     | L Re * R Re - L Im * R Im | L Re * R Im + L Im * R Re |
    //
    // It's also possible to do an fmul add sub
    // Which does (Lv1 fmul Lv2) add/sub Rv3

    // The algorithm "uses" 512 bit vectors. Since the AVX-512 instruction set
    // isn't widely available the code uses 2 256 bit vectors.

    // lhs    | L0 Re | L0 Im | L1 Re | L1 Im | L2 Re | L2 Im | L3 Re | L3 Im |
    // rhs    | R0 Re | R0 Im | R1 Re | R1 Im | R2 Re | R2 Im | R3 Re | R3 Im |

    __m256 Lv1 =
        _mm256_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b10'10'00'00);
    __m256 Rv1 =
        _mm256_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b11'11'01'01);
    __m256 Rv2 =
        _mm256_shuffle_ps(rhs.data_.Value(), rhs.data_.Value(), 0b10'11'00'01);
    __m256 Rv3 = _mm256_mul_ps(Rv1, Rv2);
    return VectorFloat8{_mm256_fmaddsub_ps(Lv1, rhs.data_.Value(), Rv3)};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorComplexFloat4
  operator*(VectorComplexFloat4 lhs, std::complex<float> rhs) noexcept {
    return lhs * VectorComplexFloat4{rhs, rhs, rhs, rhs};
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorComplexFloat4
  operator*(std::complex<float> lhs, VectorComplexFloat4 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  [[gnu::target("avx2,fma")]] friend std::ostream& operator<<(
      std::ostream& output, VectorComplexFloat4 value) {
    output << '[' << value[0] << ", " << value[1] << ", " << value[2] << ", "
           << value[3] << ']';
    return output;
  }

 private:
  VectorFloat8 data_;
};

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H
