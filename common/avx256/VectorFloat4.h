// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_FLOAT_4_H
#define AOCOMMON_AVX256_VECTOR_FLOAT_4_H

/**
 * Implements a Vector with 4 float values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexFloat2, but can be used by itself.
 *
 * @warning All functions in this header need to use a target attribute
 * like @c [[gnu::target("avx2,fma")]]. When this is not done the GCC
 * doesn't adhere to the proper ABI leading to broken code.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @note since 4 floats use 128-bits this vector uses AVX-128 vectors, despite
 * the 256 name.
 *
 * @todo Move this to aocommon when the class has matured further.
 */

#include <cassert>
#include <immintrin.h>
#include <memory>

namespace aocommon::Avx256 {

class VectorFloat4 {
 public:
  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat4() noexcept
      : data_{_mm_setzero_ps()} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] /* implicit */ VectorFloat4(
      __m128 data) noexcept
      : data_{data} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat4(
      float value) noexcept
      : data_{_mm_set1_ps(value)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat4(
      float a, float b, float c, float d) noexcept
      : data_{_mm_setr_ps(a, b, c, d)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat4(
      const float vector[4]) noexcept
      : data_{_mm_loadu_ps(vector)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] float operator[](
      size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] __m128 Value() const noexcept {
    return data_;
  }

  [[gnu::target("avx2,fma")]] VectorFloat4& operator+=(
      VectorFloat4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  [[gnu::target("avx2,fma")]] VectorFloat4& operator-=(
      VectorFloat4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat4 operator+(
      VectorFloat4 lhs, VectorFloat4 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat4 operator-(
      VectorFloat4 lhs, VectorFloat4 rhs) noexcept {
    return lhs -= rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat4 operator*(
      VectorFloat4 lhs, VectorFloat4 rhs) noexcept {
    return _mm_mul_ps(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat4 operator^(
      VectorFloat4 lhs, VectorFloat4 rhs) noexcept {
    return _mm_xor_ps(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      VectorFloat4 lhs, VectorFloat4 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3];
  }

 private:
  __m128 data_;
};

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_VECTOR_FLOAT_4_H
