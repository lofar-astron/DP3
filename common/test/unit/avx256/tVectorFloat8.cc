// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorFloat8.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(VectorFloat8)

static_assert(!std::is_default_constructible_v<aocommon::Avx256::VectorFloat8>);
static_assert(std::is_nothrow_destructible_v<aocommon::Avx256::VectorFloat8>);
static_assert(
    std::is_nothrow_copy_constructible_v<aocommon::Avx256::VectorFloat8>);
static_assert(
    std::is_nothrow_move_constructible_v<aocommon::Avx256::VectorFloat8>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorFloat8>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorFloat8>);

BOOST_AUTO_TEST_CASE(constructor_m_256) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat8, __m256>);
  const __m256 input = _mm256_set1_ps(1.0);
  const aocommon::Avx256::VectorFloat8 result{input};

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 1.0);
  BOOST_TEST(result[3] == 1.0);
  BOOST_TEST(result[4] == 1.0);
  BOOST_TEST(result[5] == 1.0);
  BOOST_TEST(result[6] == 1.0);
  BOOST_TEST(result[7] == 1.0);
}

BOOST_AUTO_TEST_CASE(constructor_float) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat8, float>);
  const aocommon::Avx256::VectorFloat8 result{-1.0};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
  BOOST_TEST(result[4] == -1.0);
  BOOST_TEST(result[5] == -1.0);
  BOOST_TEST(result[6] == -1.0);
  BOOST_TEST(result[7] == -1.0);
}

BOOST_AUTO_TEST_CASE(constructor_4_float) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat8,
                                                float, float, float, float,
                                                float, float, float, float>);
  const aocommon::Avx256::VectorFloat8 result{-1.0, 1.0,   3.75, -3.75,
                                              99.0, -99.0, 1.5,  -1.5};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == -3.75);
  BOOST_TEST(result[4] == 99.0);
  BOOST_TEST(result[5] == -99.0);
  BOOST_TEST(result[6] == 1.5);
  BOOST_TEST(result[7] == -1.5);
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat8,
                                                const float[8]>);
  const float input[] = {-1.0, 1.0, 3.75, -3.75, 99.0, -99.0, 1.5, -1.5};
  const aocommon::Avx256::VectorFloat8 result{input};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == -3.75);
  BOOST_TEST(result[4] == 99.0);
  BOOST_TEST(result[5] == -99.0);
  BOOST_TEST(result[6] == 1.5);
  BOOST_TEST(result[7] == -1.5);
}

BOOST_AUTO_TEST_CASE(constructor_double_pointer) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorFloat8,
                                                const double[8]>);
  const double input[] = {-1.0, 1.0, 3.75, -3.75, 99.0, -99.0, 1.5, -1.5};
  const aocommon::Avx256::VectorFloat8 result{input};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == -3.75);
  BOOST_TEST(result[4] == 99.0);
  BOOST_TEST(result[5] == -99.0);
  BOOST_TEST(result[6] == 1.5);
  BOOST_TEST(result[7] == -1.5);
}

BOOST_AUTO_TEST_CASE(value) {
  const aocommon::Avx256::VectorFloat8 input{-1.0};
  const aocommon::Avx256::VectorFloat8 result{input.Value()};
  static_assert(noexcept(input.Value()));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
  BOOST_TEST(result[4] == -1.0);
  BOOST_TEST(result[5] == -1.0);
  BOOST_TEST(result[6] == -1.0);
  BOOST_TEST(result[7] == -1.0);
}

BOOST_AUTO_TEST_CASE(multiply) {
  const float input[] = {-1.0, 1.0, 3.75, -3.75, 99.0, -99.0, 1.5, -1.5};
  const aocommon::Avx256::VectorFloat8 lhs{input};
  const aocommon::Avx256::VectorFloat8 rhs{input};
  const aocommon::Avx256::VectorFloat8 result = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 14.0625);
  BOOST_TEST(result[3] == 14.0625);
  BOOST_TEST(result[4] == 9801.0);
  BOOST_TEST(result[5] == 9801.0);
  BOOST_TEST(result[6] == 2.25);
  BOOST_TEST(result[7] == 2.25);
}

BOOST_AUTO_TEST_CASE(bitwise_xor) {
  const aocommon::Avx256::VectorFloat8 lhs{-1.0, 1.0,   3.75, -3.75,
                                           99.0, -99.0, 1.5,  -1.5};
  // The mask used to determine the conjugate.
  const aocommon::Avx256::VectorFloat8 rhs{0.0, -0.0, 0.0, -0.0,
                                           0.0, -0.0, 0.0, -0.0};
  const aocommon::Avx256::VectorFloat8 result = lhs ^ rhs;
  static_assert(noexcept(lhs ^ rhs));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == 3.75);
  BOOST_TEST(result[3] == 3.75);
  BOOST_TEST(result[4] == 99.0);
  BOOST_TEST(result[5] == 99.0);
  BOOST_TEST(result[6] == 1.5);
  BOOST_TEST(result[7] == 1.5);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(noexcept(aocommon::Avx256::VectorFloat8{0.0} ==
                         aocommon::Avx256::VectorFloat8{0.0}));

  BOOST_TEST(
      (aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} ==
       aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 42.0, 0.0, 0.0, 0.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 42.0, 0.0, 0.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));

  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 42.0, 0.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 42.0, 0.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 42.0, 0.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 42.0} ==
      aocommon::Avx256::VectorFloat8{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // defined(__AVX2__)
