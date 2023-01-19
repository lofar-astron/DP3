// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorComplexFloat2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

BOOST_AUTO_TEST_SUITE(VectorComplexFloat2)

static_assert(
    std::is_default_constructible_v<aocommon::Avx256::VectorComplexFloat2>);
static_assert(
    std::is_nothrow_destructible_v<aocommon::Avx256::VectorComplexFloat2>);
static_assert(std::is_nothrow_copy_constructible_v<
              aocommon::Avx256::VectorComplexFloat2>);
static_assert(std::is_nothrow_move_constructible_v<
              aocommon::Avx256::VectorComplexFloat2>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorComplexFloat2>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorComplexFloat2>);

BOOST_AUTO_TEST_CASE(constructor_default) {
  static_assert(std::is_nothrow_default_constructible_v<
                aocommon::Avx256::VectorComplexFloat2>);
  const aocommon::Avx256::VectorComplexFloat2 result;

  BOOST_TEST(result[0] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(constructor_vector_float) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat2,
                                      aocommon::Avx256::VectorFloat4>);
  const __m128 input = _mm_set1_ps(1.0);
  const aocommon::Avx256::VectorComplexFloat2 result{
      aocommon::Avx256::VectorFloat4{input}};

  BOOST_TEST(result[0] == (std::complex<float>{1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{1.0, 1.0}));
}

BOOST_AUTO_TEST_CASE(constructor_2_complex_float) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat2,
                                      std::complex<float>,
                                      std::complex<float>>);

  const aocommon::Avx256::VectorComplexFloat2 result{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75}};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat2,
                                      const std::complex<float>[2]>);
  const std::complex<float> input[] = {{-1.0, 1.0}, {3.75, -3.75}};
  const aocommon::Avx256::VectorComplexFloat2 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(operator_m128) {
  static_assert(
      noexcept(static_cast<__m128>(aocommon::Avx256::VectorComplexFloat2{
          static_cast<const std::complex<float>*>(nullptr)})));

  const aocommon::Avx256::VectorComplexFloat2 input{
      aocommon::Avx256::VectorFloat4{-1.0, 1.0, 3.75, -3.75}};

  BOOST_TEST(static_cast<__m128>(input)[0] == -1.0);
  BOOST_TEST(static_cast<__m128>(input)[1] == 1.0);
  BOOST_TEST(static_cast<__m128>(input)[2] == 3.75);
  BOOST_TEST(static_cast<__m128>(input)[3] == -3.75);
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(aocommon::Avx256::VectorComplexFloat2{
      static_cast<const std::complex<float>*>(nullptr)}
                             .Conjugate()));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{1.0, 2.0}, {10, 11}}.Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{1.0, -2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat2{{-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(operator_plus_minus) {
  aocommon::Avx256::VectorComplexFloat2 r{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexFloat2 value{{4, 8}, {40, 44}};

  r -= value;
  static_assert(noexcept(r -= value));

  BOOST_TEST(r[0] == (std::complex<float>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<float>{-30, -33}));
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  aocommon::Avx256::VectorComplexFloat2 r{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexFloat2 value{{4, 8}, {40, 44}};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r[0] == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<float>{50.0, 55.0}));
}

BOOST_AUTO_TEST_CASE(operator_plus) {
  const aocommon::Avx256::VectorComplexFloat2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexFloat2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexFloat2 r = lhs + rhs;
  static_assert(noexcept(lhs + rhs));

  BOOST_TEST(r[0] == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<float>{50.0, 55.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus) {
  const aocommon::Avx256::VectorComplexFloat2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexFloat2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexFloat2 r = lhs - rhs;
  static_assert(noexcept(lhs - rhs));

  BOOST_TEST(r[0] == (std::complex<float>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<float>{-30, -33}));
}

BOOST_AUTO_TEST_CASE(multiply) {
  const aocommon::Avx256::VectorComplexFloat2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexFloat2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexFloat2 r = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_CHECK_CLOSE(r[0].real(), -12.0, 1e-6);
  BOOST_CHECK_CLOSE(r[0].imag(), 16.0, 1e-6);
  BOOST_CHECK_CLOSE(r[1].real(), -84.0, 1e-6);
  BOOST_CHECK_CLOSE(r[1].imag(), 880.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(aocommon::Avx256::VectorComplexFloat2{
                   static_cast<const std::complex<float>*>(nullptr)} ==
               aocommon::Avx256::VectorComplexFloat2{
                   static_cast<const std::complex<float>*>(nullptr)}));

  BOOST_TEST((aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}} ==
              aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat2{{42.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat2{{0.0, 42.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {42.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 42.0}} ==
               aocommon::Avx256::VectorComplexFloat2{{0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_CASE(output) {
  const aocommon::Avx256::VectorComplexFloat2 input{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75}};

  std::stringstream result;
  result << input;

  BOOST_CHECK_EQUAL(result.str(), "[(-1,1), (3.75,-3.75)]");
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // defined(__AVX2__)
