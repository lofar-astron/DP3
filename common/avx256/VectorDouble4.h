// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
#define AOCOMMON_AVX256_VECTOR_DOUBLE_4_H

/**
 * @file Implements a Vector with 4 double values.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexDouble2, but can be used by itself.
 *
 * @note The class only implements a subset of the vector operations. Other
 * operations will be added on a when-needed basis.
 *
 * @todo Move this to aocommon then the class has matured further.
 */
#include <cassert>
#include <memory>

#if defined(__AVX2__)

#include <immintrin.h>

namespace aocommon::Avx256 {

class VectorDouble4 {
 public:
  /* implicit */ VectorDouble4(__m256d data) noexcept : data_{data} {}

  explicit VectorDouble4(double value) noexcept
      : data_{_mm256_set1_pd(value)} {}

  explicit VectorDouble4(double a, double b, double c, double d) noexcept
      : data_{_mm256_setr_pd(a, b, c, d)} {}

  explicit VectorDouble4(const double vector[4]) noexcept
      : data_{_mm256_loadu_pd(vector)} {}

  explicit VectorDouble4(const float vector[4]) noexcept
      : data_(_mm256_cvtps_pd(_mm_loadu_ps(std::addressof(vector[0])))) {}

  [[nodiscard]] double operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  [[nodiscard]] __m256d Value() const noexcept { return data_; }

  VectorDouble4& operator+=(VectorDouble4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  VectorDouble4& operator-=(VectorDouble4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  [[nodiscard]] friend VectorDouble4 operator+(VectorDouble4 lhs,
                                               VectorDouble4 rhs) noexcept {
    return lhs += rhs;
  }

  [[nodiscard]] friend VectorDouble4 operator-(VectorDouble4 lhs,
                                               VectorDouble4 rhs) noexcept {
    return lhs -= rhs;
  }

  [[nodiscard]] friend VectorDouble4 operator*(VectorDouble4 lhs,
                                               VectorDouble4 rhs) noexcept {
    return _mm256_mul_pd(lhs.data_, rhs.data_);
  }

  [[nodiscard]] friend VectorDouble4 operator^(VectorDouble4 lhs,
                                               VectorDouble4 rhs) noexcept {
    return _mm256_xor_pd(lhs.data_, rhs.data_);
  }

  [[nodiscard]] friend bool operator==(VectorDouble4 lhs,
                                       VectorDouble4 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3];
  }

 private:
  __m256d data_;
};

}  // namespace aocommon::Avx256

#endif  // defined(__AVX2__)

#endif  // AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
