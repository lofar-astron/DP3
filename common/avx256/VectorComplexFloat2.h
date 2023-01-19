// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_2_H
#define AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_2_H

/**
 * @file Implements a Vector with 2 complex float values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::MatrixComplexFloat2x2, but can be used by itself.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */

#include "common/avx256/VectorFloat4.h"

#include <cassert>
#include <complex>
#include <ostream>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class VectorComplexFloat2 {
 public:
  [[nodiscard]] VectorComplexFloat2() noexcept = default;

  /* implicit */ VectorComplexFloat2(VectorFloat4 data) noexcept
      : data_{data} {}

  explicit VectorComplexFloat2(std::complex<float> a,
                               std::complex<float> b) noexcept
      : data_{VectorFloat4{a.real(), a.imag(), b.real(), b.imag()}} {}

  explicit VectorComplexFloat2(const std::complex<float> vector[2]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat4{
            reinterpret_cast<const float*>(std::addressof(vector[0]))}} {}

  [[nodiscard]] std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return {data_[2 * index], data_[2 * index + 1]};
  }

  [[nodiscard]] explicit operator __m128() const noexcept {
    return data_.Value();
  }

  [[nodiscard]] VectorComplexFloat2 Conjugate() const noexcept {
    // Xor-ing a float with  0.0 will not change the value.
    // Xor-ing a float with -0.0 will change the sign of the value.
    __m128 mask = _mm_setr_ps(0.0, -0.0, 0.0, -0.0);
    return data_ ^ mask;
  }

  VectorComplexFloat2& operator+=(VectorComplexFloat2 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  VectorComplexFloat2& operator-=(VectorComplexFloat2 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] friend VectorComplexFloat2 operator+(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] friend VectorComplexFloat2 operator-(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
    return lhs -= rhs;
  }

  /// Multiplies the elements of 2 vectors on parallel.
  ///
  /// r[0] = lhs[0] * rhs[0]
  /// r[1] = lhs[1] * rhs[1]
  [[nodiscard]] friend VectorComplexFloat2 operator*(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
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

    // The algorithm "uses" 256 bit vectors in two 128 bit vectors. Testing
    // with 1 256 bit vector didn't improve much. The main issue is:
    //   128 bit: mul + fmaddsub
    //   256 bit: mul + addsub
    // The performace of fmaddsub and addsub seem to be the same according to
    // Intel. The additional register manimulation actually makes it more
    // expense.
    //
    // Not there can be an improvement by doing two diagonal matrixes in
    // parallel. However that requires more code changes.

    __m128 Rv1 =
        _mm_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b11'11'01'01);
    __m128 Rv2 =
        _mm_shuffle_ps(rhs.data_.Value(), rhs.data_.Value(), 0b10'11'00'01);
    __m128 Lv1 =
        _mm_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b10'10'00'00);
    __m128 Rv3 = _mm_mul_ps(Rv1, Rv2);
    return VectorFloat4{_mm_fmaddsub_ps(Lv1, rhs.data_.Value(), Rv3)};
  }

  friend VectorComplexFloat2 operator*(VectorComplexFloat2 lhs,
                                       std::complex<float> rhs) noexcept {
    return lhs * VectorComplexFloat2{rhs, rhs};
  }

  friend VectorComplexFloat2 operator*(std::complex<float> lhs,
                                       VectorComplexFloat2 rhs) noexcept {
    return rhs * lhs;
  }

  [[nodiscard]] friend bool operator==(VectorComplexFloat2 lhs,
                                       VectorComplexFloat2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  friend std::ostream& operator<<(std::ostream& output,
                                  VectorComplexFloat2 value) {
    output << '[' << value[0] << ", " << value[1] << ']';
    return output;
  }

 private:
  VectorFloat4 data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_2_H
