// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_FLOAT_8_H
#define AOCOMMON_AVX256_VECTOR_FLOAT_8_H

/**
 * Implements a Vector with 8 float values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexFloat4, but can be used by itself.
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

#include <cassert>
#include <immintrin.h>
#include <memory>

namespace aocommon::Avx256 {

class VectorFloat8 {
 public:
  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat8() noexcept
      : data_{_mm256_setzero_ps()} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] /* implicit */ VectorFloat8(
      __m256 data) noexcept
      : data_{data} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat8(
      float value) noexcept
      : data_{_mm256_set1_ps(value)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat8(
      float a, float b, float c, float d, float e, float f, float g,
      float h) noexcept
      : data_{_mm256_setr_ps(a, b, c, d, e, f, g, h)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat8(
      const float vector[8]) noexcept
      : data_{_mm256_loadu_ps(vector)} {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] explicit VectorFloat8(
      const double vector[8]) noexcept {
    __m256 lo = _mm256_castps128_ps256(
        _mm256_cvtpd_ps(_mm256_loadu_pd(std::addressof(vector[0]))));
    __m128 hi = _mm256_cvtpd_ps(_mm256_loadu_pd(std::addressof(vector[4])));

    data_ = _mm256_insertf128_ps(lo, hi, 1);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] float operator[](
      size_t index) const noexcept {
    assert(index < 8 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] __m256 Value() const noexcept {
    return data_;
  }

  /** Assign data stored by 8 element vector to destination buffer */
  [[gnu::target("avx2,fma")]] void AssignTo(float* destination) const noexcept {
    _mm256_storeu_ps(destination, data_);
  }

  [[gnu::target("avx2,fma")]] VectorFloat8& operator+=(
      VectorFloat8 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  [[gnu::target("avx2,fma")]] VectorFloat8& operator-=(
      VectorFloat8 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat8 operator+(
      VectorFloat8 lhs, VectorFloat8 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat8 operator-(
      VectorFloat8 lhs, VectorFloat8 rhs) noexcept {
    return lhs -= rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat8 operator*(
      VectorFloat8 lhs, VectorFloat8 rhs) noexcept {
    return _mm256_mul_ps(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorFloat8 operator^(
      VectorFloat8 lhs, VectorFloat8 rhs) noexcept {
    return _mm256_xor_ps(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      VectorFloat8 lhs, VectorFloat8 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3] &&
           lhs.data_[4] == rhs.data_[4] && lhs.data_[5] == rhs.data_[5] &&
           lhs.data_[6] == rhs.data_[6] && lhs.data_[7] == rhs.data_[7];
  }

 private:
  __m256 data_;
};

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_VECTOR_FLOAT_8_H
