// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorDouble4.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// Silences the diagnostic "ignoring attributes on template argument ‘__m256d’"
// This diagnostic is issued on code that's part of libstdc++.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(VectorDouble4)

static_assert(
    !std::is_default_constructible_v<aocommon::Avx256::VectorDouble4>);
static_assert(std::is_nothrow_destructible_v<aocommon::Avx256::VectorDouble4>);
static_assert(
    std::is_nothrow_copy_constructible_v<aocommon::Avx256::VectorDouble4>);
static_assert(
    std::is_nothrow_move_constructible_v<aocommon::Avx256::VectorDouble4>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorDouble4>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorDouble4>);

BOOST_AUTO_TEST_CASE(constructor_m_256) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorDouble4,
                                                __m256d>);
  const __m256d input = _mm256_set1_pd(1.0);
  const aocommon::Avx256::VectorDouble4 result{input};

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 1.0);
  BOOST_TEST(result[3] == 1.0);
}

BOOST_AUTO_TEST_CASE(constructor_double) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorDouble4, double>);
  const aocommon::Avx256::VectorDouble4 result{-1.0};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
}

BOOST_AUTO_TEST_CASE(constructor_4_double) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorDouble4, double,
                                      double, double, double>);
  const aocommon::Avx256::VectorDouble4 result{-1.0, 1.0, 4.2, -4.2};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 4.2);
  BOOST_TEST(result[3] == -4.2);
}

BOOST_AUTO_TEST_CASE(constructor_double_pointer) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorDouble4,
                                                const double[4]>);
  const double input[] = {-1.0, 1.0, 4.2, -4.2};
  const aocommon::Avx256::VectorDouble4 result{input};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 4.2);
  BOOST_TEST(result[3] == -4.2);
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(std::is_nothrow_constructible_v<aocommon::Avx256::VectorDouble4,
                                                const float[4]>);
  const float input[] = {-1.0, 1.0, 2.5, -2.5};
  const aocommon::Avx256::VectorDouble4 result{input};

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 2.5);
  BOOST_TEST(result[3] == -2.5);
}

BOOST_AUTO_TEST_CASE(value) {
  const aocommon::Avx256::VectorDouble4 input{-1.0};
  const aocommon::Avx256::VectorDouble4 result{input.Value()};
  static_assert(noexcept(input.Value()));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == -1.0);
  BOOST_TEST(result[3] == -1.0);
}

BOOST_AUTO_TEST_CASE(multiply) {
  const double input[] = {-1.0, 1.0, 4.2, -4.2};
  const aocommon::Avx256::VectorDouble4 lhs{input};
  const aocommon::Avx256::VectorDouble4 rhs{input};
  const aocommon::Avx256::VectorDouble4 result = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(result[0] == 1.0);
  BOOST_TEST(result[1] == 1.0);
  BOOST_TEST(result[2] == 17.64);
  BOOST_TEST(result[3] == 17.64);
}

BOOST_AUTO_TEST_CASE(bitwise_xor) {
  const aocommon::Avx256::VectorDouble4 lhs{-1.0, 1.0, 4.2, -4.2};
  // The mask used to determine the conjugate.
  const aocommon::Avx256::VectorDouble4 rhs{0.0, -0.0, 0.0, -0.0};
  const aocommon::Avx256::VectorDouble4 result = lhs ^ rhs;
  static_assert(noexcept(lhs ^ rhs));

  BOOST_TEST(result[0] == -1.0);
  BOOST_TEST(result[1] == -1.0);
  BOOST_TEST(result[2] == 4.2);
  BOOST_TEST(result[3] == 4.2);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(noexcept(aocommon::Avx256::VectorDouble4{0.0} ==
                         aocommon::Avx256::VectorDouble4{0.0}));

  BOOST_TEST((aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0} ==
              aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorDouble4{42.0, 0.0, 0.0, 0.0} ==
               aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorDouble4{0.0, 42.0, 0.0, 0.0} ==
               aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorDouble4{0.0, 0.0, 42.0, 0.0} ==
               aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0}));
  BOOST_TEST(!(aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 42.0} ==
               aocommon::Avx256::VectorDouble4{0.0, 0.0, 0.0, 0.0}));
}

BOOST_AUTO_TEST_SUITE_END()

#pragma GCC diagnostic pop

#endif  // defined(__AVX2__)
