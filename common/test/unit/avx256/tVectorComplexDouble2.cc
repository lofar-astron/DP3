// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorComplexDouble2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(VectorComplexDouble2)

static_assert(
    std::is_default_constructible_v<aocommon::Avx256::VectorComplexDouble2>);
static_assert(
    std::is_nothrow_destructible_v<aocommon::Avx256::VectorComplexDouble2>);
static_assert(std::is_nothrow_copy_constructible_v<
              aocommon::Avx256::VectorComplexDouble2>);
static_assert(std::is_nothrow_move_constructible_v<
              aocommon::Avx256::VectorComplexDouble2>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorComplexDouble2>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorComplexDouble2>);

BOOST_AUTO_TEST_CASE(constructor_default) {
  static_assert(std::is_nothrow_default_constructible_v<
                aocommon::Avx256::VectorComplexDouble2>);
  const aocommon::Avx256::VectorComplexDouble2 result;

  BOOST_TEST(result[0] == (std::complex<double>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<double>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(constructor_vector_double) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexDouble2,
                                      aocommon::Avx256::VectorDouble4>);
  const __m256d input = _mm256_set1_pd(1.0);
  const aocommon::Avx256::VectorComplexDouble2 result{
      aocommon::Avx256::VectorDouble4{input}};

  BOOST_TEST(result[0] == (std::complex<double>{1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<double>{1.0, 1.0}));
}

BOOST_AUTO_TEST_CASE(constructor_2_complex_double) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexDouble2,
                                      std::complex<double>,
                                      std::complex<double>>);

  const aocommon::Avx256::VectorComplexDouble2 result{
      std::complex<double>{-1.0, 1.0}, std::complex<double>{3.75, -3.75}};

  BOOST_TEST(result[0] == (std::complex<double>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<double>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexDouble2,
                                      const std::complex<float>[2]>);
  const std::complex<float> input[] = {{-1.0, 1.0}, {3.75, -3.75}};
  const aocommon::Avx256::VectorComplexDouble2 result{input};

  BOOST_TEST(result[0] == (std::complex<double>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<double>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(constructor_double_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexDouble2,
                                      const std::complex<double>[2]>);
  const std::complex<double> input[] = {{-1.0, 1.0}, {3.75, -3.75}};
  const aocommon::Avx256::VectorComplexDouble2 result{input};

  BOOST_TEST(result[0] == (std::complex<double>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<double>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(operator_m256d) {
  static_assert(
      noexcept(static_cast<__m256d>(aocommon::Avx256::VectorComplexDouble2{
          static_cast<const std::complex<double>*>(nullptr)})));

  const aocommon::Avx256::VectorComplexDouble2 input{
      aocommon::Avx256::VectorDouble4{-1.0, 1.0, 3.75, -3.75}};

  BOOST_TEST(static_cast<__m256d>(input)[0] == -1.0);
  BOOST_TEST(static_cast<__m256d>(input)[1] == 1.0);
  BOOST_TEST(static_cast<__m256d>(input)[2] == 3.75);
  BOOST_TEST(static_cast<__m256d>(input)[3] == -3.75);
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(aocommon::Avx256::VectorComplexDouble2{
      static_cast<const std::complex<double>*>(nullptr)}
                             .Conjugate()));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{1.0, 2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{1.0, -2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, -2.0}, {-10, -11}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexDouble2{{-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(operator_plus_minus) {
  aocommon::Avx256::VectorComplexDouble2 r{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 value{{4, 8}, {40, 44}};

  r -= value;
  static_assert(noexcept(r -= value));

  BOOST_TEST(r[0] == (std::complex<double>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<double>{-30, -33}));
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  aocommon::Avx256::VectorComplexDouble2 r{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 value{{4, 8}, {40, 44}};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r[0] == (std::complex<double>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<double>{50.0, 55.0}));
}

BOOST_AUTO_TEST_CASE(operator_plus) {
  const aocommon::Avx256::VectorComplexDouble2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexDouble2 r = lhs + rhs;
  static_assert(noexcept(lhs + rhs));

  BOOST_TEST(r[0] == (std::complex<double>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<double>{50.0, 55.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus) {
  const aocommon::Avx256::VectorComplexDouble2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexDouble2 r = lhs - rhs;
  static_assert(noexcept(lhs - rhs));

  BOOST_TEST(r[0] == (std::complex<double>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<double>{-30, -33}));
}

BOOST_AUTO_TEST_CASE(multiply) {
  const aocommon::Avx256::VectorComplexDouble2 lhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 rhs{{4, 8}, {40, 44}};

  const aocommon::Avx256::VectorComplexDouble2 r = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(r[0] == (std::complex<double>{-12.0, 16.0}));
  BOOST_TEST(r[1] == (std::complex<double>{-84.0, 880.0}));
}

BOOST_AUTO_TEST_CASE(multiply_vector_and_value) {
  const aocommon::Avx256::VectorComplexDouble2 lhs{{1.0, 2.0}, {10, 11}};

  const std::complex<double> rhs{4, 8};

  const aocommon::Avx256::VectorComplexDouble2 r = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(r[0] == (std::complex<double>{-12.0, 16.0}));
  BOOST_TEST(r[1] == (std::complex<double>{-48, 124}));
}

BOOST_AUTO_TEST_CASE(multiply_value_and_vector) {
  const std::complex<double> lhs{4, 8};

  const aocommon::Avx256::VectorComplexDouble2 rhs{{1.0, 2.0}, {10, 11}};

  const aocommon::Avx256::VectorComplexDouble2 r = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(r[0] == (std::complex<double>{-12.0, 16.0}));
  BOOST_TEST(r[1] == (std::complex<double>{-48, 124}));
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(aocommon::Avx256::VectorComplexDouble2{
                   static_cast<const std::complex<double>*>(nullptr)} ==
               aocommon::Avx256::VectorComplexDouble2{
                   static_cast<const std::complex<double>*>(nullptr)}));

  BOOST_TEST((aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}} ==
              aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(
      !(aocommon::Avx256::VectorComplexDouble2{{42.0, 0.0}, {0.0, 0.0}} ==
        aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(
      !(aocommon::Avx256::VectorComplexDouble2{{0.0, 42.0}, {0.0, 0.0}} ==
        aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(
      !(aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {42.0, 0.0}} ==
        aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(
      !(aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 42.0}} ==
        aocommon::Avx256::VectorComplexDouble2{{0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_CASE(output) {
  const aocommon::Avx256::VectorComplexDouble2 input{
      std::complex<double>{-1.0, 1.0}, std::complex<double>{3.75, -3.75}};

  std::stringstream result;
  result << input;

  BOOST_CHECK_EQUAL(result.str(), "[(-1,1), (3.75,-3.75)]");
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // defined(__AVX2__)
