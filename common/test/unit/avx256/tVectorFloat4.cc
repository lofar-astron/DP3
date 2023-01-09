// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorFloat4.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// Silences the diagnostic "ignoring attributes on template argument ‘_mm128’"
// This diagnostic is issued on code that's part of libstdc++.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(VectorFloat4)

static_assert(std::is_default_constructible_v<aocommon::Avx256::VectorFloat4>);
static_assert(std::is_nothrow_destructible_v<aocommon::Avx256::VectorFloat4>);
static_assert(
    std::is_nothrow_copy_constructible_v<aocommon::Avx256::VectorFloat4>);
static_assert(
    std::is_nothrow_move_constructible_v<aocommon::Avx256::VectorFloat4>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorFloat4>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorFloat4>);

BOOST_AUTO_TEST_CASE(constructor_default) {
  static_assert(
      std::is_nothrow_default_constructible_v<aocommon::Avx256::VectorFloat4>);
  const aocommon::Avx256::VectorFloat4 result;

  BOOST_TEST(result[0] == 0.0);
  BOOST_TEST(result[1] == 0.0);
  BOOST_TEST(result[2] == 0.0);
  BOOST_TEST(result[3] == 0.0);
}

BOOST_AUTO_TEST_CASE(constructor_m_128) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat4, __m128>);
  const __m128 input = _mm_set1_ps(1.0);
  const aocommon::Avx256::VectorFloat4 result{input};

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 1.0);
  BOOST_TEST(result[3] == 1.0);
}

BOOST_AUTO_TEST_CASE(constructor_float) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat4, float>);
  const aocommon::Avx256::VectorFloat4 result{-1.0};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
}

BOOST_AUTO_TEST_CASE(constructor_4_float) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat4,
                                                float, float, float, float>);
  const aocommon::Avx256::VectorFloat4 result{-1.0, 1.0, 3.75, -3.75};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == -3.75);
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat4,
                                                const float[4]>);
  const float input[] = {-1.0, 1.0, 3.75, -3.75};
  const aocommon::Avx256::VectorFloat4 result{input};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == -3.75);
}

BOOST_AUTO_TEST_CASE(value) {
  const aocommon::Avx256::VectorFloat4 input{-1.0};
  const aocommon::Avx256::VectorFloat4 result{input.Value()};
  static_assert(noexcept(input.Value()));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  aocommon::Avx256::VectorFloat4 r{1, 2, 10, 11};

  const aocommon::Avx256::VectorFloat4 value{4, 8, 40, 44};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r[0] == 5);
  BOOST_TEST(r[1] == 10);
  BOOST_TEST(r[2] == 50);
  BOOST_TEST(r[3] == 55);
}

BOOST_AUTO_TEST_CASE(operator_minus_equal) {
  aocommon::Avx256::VectorFloat4 r{1, 2, 10, 11};

  const aocommon::Avx256::VectorFloat4 value{4, 8, 40, 44};

  r -= value;
  static_assert(noexcept(r -= value));

  BOOST_TEST(r[0] == -3);
  BOOST_TEST(r[1] == -6);
  BOOST_TEST(r[2] == -30);
  BOOST_TEST(r[3] == -33);
}

BOOST_AUTO_TEST_CASE(operator_plus) {
  const aocommon::Avx256::VectorFloat4 lhs{1.0, 2.0, 10, 11};

  const aocommon::Avx256::VectorFloat4 rhs{4, 8, 40, 44};

  const aocommon::Avx256::VectorFloat4 r = lhs + rhs;
  static_assert(noexcept(lhs + rhs));

  BOOST_TEST(r[0] == 5);
  BOOST_TEST(r[1] == 10);
  BOOST_TEST(r[2] == 50);
  BOOST_TEST(r[3] == 55);
}

BOOST_AUTO_TEST_CASE(operator_minus) {
  const aocommon::Avx256::VectorFloat4 lhs{1.0, 2.0, 10, 11};

  const aocommon::Avx256::VectorFloat4 rhs{4, 8, 40, 44};

  const aocommon::Avx256::VectorFloat4 r = lhs - rhs;
  static_assert(noexcept(lhs - rhs));

  BOOST_TEST(r[0] == -3);
  BOOST_TEST(r[1] == -6);
  BOOST_TEST(r[2] == -30);
  BOOST_TEST(r[3] == -33);
}

BOOST_AUTO_TEST_CASE(multiply) {
  const float input[] = {-1.0, 1.0, 3.75, -3.75};
  const aocommon::Avx256::VectorFloat4 lhs{input};
  const aocommon::Avx256::VectorFloat4 rhs{input};
  const aocommon::Avx256::VectorFloat4 result = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 14.0625);
  BOOST_TEST(result[3] == 14.0625);
}

BOOST_AUTO_TEST_CASE(bitwise_xor) {
  const aocommon::Avx256::VectorFloat4 lhs{-1.0, 1.0, 3.75, -3.75};
  // The mask used to determine the conjugate.
  const aocommon::Avx256::VectorFloat4 rhs{0.0, -0.0, 0.0, -0.0};
  const aocommon::Avx256::VectorFloat4 result = lhs ^ rhs;
  static_assert(noexcept(lhs ^ rhs));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == 3.75);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(noexcept(aocommon::Avx256::VectorFloat4{0.0} ==
                         aocommon::Avx256::VectorFloat4{0.0}));

  BOOST_TEST((aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0} ==
              aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorFloat4{42.0, 0.0, 0.0, 0.0} ==
               aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorFloat4{0.0, 42.0, 0.0, 0.0} ==
               aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorFloat4{0.0, 0.0, 42.0, 0.0} ==
               aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 42.0} ==
               aocommon::Avx256::VectorFloat4{0.0, 0.0, 0.0, 0.0}));
}

BOOST_AUTO_TEST_SUITE_END()

#pragma GCC diagnostic pop

#endif  // defined(__AVX2__)
