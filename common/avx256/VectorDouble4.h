// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
#define AOCOMMON_AVX256_VECTOR_DOUBLE_4_H

/**
 * Implements a Vector with 4 double values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexDouble2, but can be used by itself.
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

class VectorDouble4 {
 public:
  [[gnu::target("avx2,fma")]] explicit VectorDouble4() noexcept
      : data_{_mm256_setzero_pd()} {}

  [[gnu::target("avx2,fma")]] [[gnu::target(
      "avx2,fma")]] /* implicit */ VectorDouble4(__m256d data) noexcept
      : data_{data} {}

  [[gnu::target("avx2,fma")]] explicit VectorDouble4(double value) noexcept
      : data_{_mm256_set1_pd(value)} {}

  [[gnu::target("avx2,fma")]] explicit VectorDouble4(double a, double b,
                                                     double c,
                                                     double d) noexcept
      : data_{_mm256_setr_pd(a, b, c, d)} {}

  [[gnu::target("avx2,fma")]] explicit VectorDouble4(
      const double vector[4]) noexcept
      : data_{_mm256_loadu_pd(vector)} {}

  [[gnu::target("avx2,fma")]] explicit VectorDouble4(
      const float vector[4]) noexcept
      : data_(_mm256_cvtps_pd(_mm_loadu_ps(std::addressof(vector[0])))) {}

  [[nodiscard]] [[gnu::target("avx2,fma")]] double operator[](
      size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] __m256d Value() const noexcept {
    return data_;
  }

  /** Assign data stored by 4 element vector to destination buffer. */
  [[gnu::target("avx2,fma")]] void AssignTo(
      double* destination) const noexcept {
    _mm256_storeu_pd(destination, data_);
  }

  [[gnu::target("avx2,fma")]] VectorDouble4& operator+=(
      VectorDouble4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  [[gnu::target("avx2,fma")]] VectorDouble4& operator-=(
      VectorDouble4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorDouble4 operator+(
      VectorDouble4 lhs, VectorDouble4 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorDouble4 operator-(
      VectorDouble4 lhs, VectorDouble4 rhs) noexcept {
    return lhs -= rhs;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorDouble4 operator*(
      VectorDouble4 lhs, VectorDouble4 rhs) noexcept {
    return _mm256_mul_pd(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend VectorDouble4 operator^(
      VectorDouble4 lhs, VectorDouble4 rhs) noexcept {
    return _mm256_xor_pd(lhs.data_, rhs.data_);
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] friend bool operator==(
      VectorDouble4 lhs, VectorDouble4 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3];
  }

 private:
  __m256d data_;
};

}  // namespace aocommon::Avx256

#endif  // AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
