// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_COMPLEX_DOUBLE_2_H
#define AOCOMMON_AVX256_VECTOR_COMPLEX_DOUBLE_2_H

/**
 * @file Implements a Vector with 2 complex double values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::MatrixComplexDouble2x2, but can be used by itself.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */

#include "common/avx256/VectorDouble4.h"

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class VectorComplexDouble2 {
 public:
  /* implicit */ VectorComplexDouble2(VectorDouble4 data) noexcept
      : data_{data} {}

  explicit VectorComplexDouble2(std::complex<double> a,
                                std::complex<double> b) noexcept
      : data_{VectorDouble4{a.real(), a.imag(), b.real(), b.imag()}} {}

  explicit VectorComplexDouble2(const std::complex<float> vector[2]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorDouble4{
            reinterpret_cast<const float*>(std::addressof(vector[0]))}} {}

  explicit VectorComplexDouble2(const std::complex<double> vector[2]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorDouble4{
            reinterpret_cast<const double*>(std::addressof(vector[0]))}} {}

  [[nodiscard]] std::complex<double> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return {data_[2 * index], data_[2 * index + 1]};
  }

  [[nodiscard]] explicit operator __m256d() const noexcept {
    return data_.Value();
  }

  [[nodiscard]] VectorComplexDouble2 Conjugate() const noexcept {
    // Xor-ing a double with  0.0 will not change the value.
    // Xor-ing a double with -0.0 will change the sign of the value.
    __m256d mask = _mm256_setr_pd(0.0, -0.0, 0.0, -0.0);
    return data_ ^ mask;
  }

  VectorComplexDouble2& operator+=(VectorComplexDouble2 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  VectorComplexDouble2& operator-=(VectorComplexDouble2 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] friend VectorComplexDouble2 operator+(
      VectorComplexDouble2 lhs, VectorComplexDouble2 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] friend VectorComplexDouble2 operator-(
      VectorComplexDouble2 lhs, VectorComplexDouble2 rhs) noexcept {
    return lhs -= rhs;
  }

  /// Multiplies the elements of 2 vectors on parallel.
  ///
  /// r[0] = lhs[0] * rhs[0]
  /// r[1] = lhs[1] * rhs[1]
  [[nodiscard]] friend VectorComplexDouble2 operator*(
      VectorComplexDouble2 lhs, VectorComplexDouble2 rhs) noexcept {
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

    // lhs    | L0 Re | L0 Im | L1 Re | L1 Im |

    __m256d Lv1 =
        _mm256_shuffle_pd(lhs.data_.Value(), lhs.data_.Value(), 0b00'00);
    __m256d Rv1 =
        _mm256_shuffle_pd(lhs.data_.Value(), lhs.data_.Value(), 0b11'11);
    __m256d Rv2 =
        _mm256_shuffle_pd(rhs.data_.Value(), rhs.data_.Value(), 0b01'01);
    __m256d Rv3 = _mm256_mul_pd(Rv1, Rv2);
    return VectorDouble4{_mm256_fmaddsub_pd(Lv1, rhs.data_.Value(), Rv3)};
  }

  [[nodiscard]] friend bool operator==(VectorComplexDouble2 lhs,
                                       VectorComplexDouble2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  VectorComplexDouble2 value) {
    output << '[' << value[0] << ", " << value[1] << ']';
    return output;
  }

 private:
  VectorDouble4 data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_VECTOR_COMPLEX_DOUBLE_2_H
