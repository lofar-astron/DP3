// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H
#define AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H

/**
 * @file Implements a Vector with 4 complex float values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::MaxtrixComplexFloat2x2, but can be used by itself.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */

#include "common/avx256/VectorFloat8.h"

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class VectorComplexFloat4 {
 public:
  /* implicit */ VectorComplexFloat4(VectorFloat8 data) noexcept
      : data_{data} {}

  explicit VectorComplexFloat4(std::complex<float> a, std::complex<float> b,
                               std::complex<float> c,
                               std::complex<float> d) noexcept
      : data_{VectorFloat8{a.real(), a.imag(), b.real(), b.imag(), c.real(),
                           c.imag(), d.real(), d.imag()}} {}

  explicit VectorComplexFloat4(const std::complex<float> vector[4]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat8{
            reinterpret_cast<const float*>(std::addressof(vector[0]))}} {}

  explicit VectorComplexFloat4(const std::complex<double> vector[4]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat8{
            reinterpret_cast<const double*>(std::addressof(vector[0]))}} {}

  [[nodiscard]] std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return {data_[2 * index], data_[2 * index + 1]};
  }

  [[nodiscard]] VectorComplexFloat4 Conjugate() const noexcept {
    // Xor-ing a float with  0.0 will not change the value.
    // Xor-ing a float with -0.0 will change the sign of the value.
    __m256 mask = _mm256_setr_ps(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0);
    return data_ ^ mask;
  }

  /// Multiplies the elements of 2 vectors on parallel.
  ///
  /// r[0] = lhs[0] * rhs[0]
  /// r[1] = lhs[1] * rhs[1]
  /// r[2] = lhs[2] * rhs[2]
  /// r[3] = lhs[3] * rhs[3]
  [[nodiscard]] friend VectorComplexFloat4 operator*(
      VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
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

  [[nodiscard]] friend VectorComplexFloat4 operator+(
      VectorComplexFloat4 lhs, VectorComplexFloat4 rhs) noexcept {
    return VectorFloat8{_mm256_add_ps(lhs.data_.Value(), rhs.data_.Value())};
  }

  [[nodiscard]] friend bool operator==(VectorComplexFloat4 lhs,
                                       VectorComplexFloat4 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  VectorComplexFloat4 value) {
    output << '[' << value[0] << ", " << value[1] << ", " << value[2] << ", "
           << value[3] << ']';
    return output;
  }

 private:
  VectorFloat8 data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_4_H
