// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorComplexFloat4.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(VectorComplexFloat4)

static_assert(
    std::is_default_constructible_v<aocommon::Avx256::VectorComplexFloat4>);
static_assert(
    std::is_nothrow_destructible_v<aocommon::Avx256::VectorComplexFloat4>);
static_assert(std::is_nothrow_copy_constructible_v<
              aocommon::Avx256::VectorComplexFloat4>);
static_assert(std::is_nothrow_move_constructible_v<
              aocommon::Avx256::VectorComplexFloat4>);
static_assert(
    std::is_nothrow_copy_assignable_v<aocommon::Avx256::VectorComplexFloat4>);
static_assert(
    std::is_nothrow_move_assignable_v<aocommon::Avx256::VectorComplexFloat4>);

BOOST_AUTO_TEST_CASE(constructor_default) {
  static_assert(std::is_nothrow_default_constructible_v<
                aocommon::Avx256::VectorComplexFloat4>);
  const aocommon::Avx256::VectorComplexFloat4 result;

  BOOST_TEST(result[0] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[2] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[3] == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(constructor_vector_float) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat4,
                                      aocommon::Avx256::VectorFloat8>);
  const __m256 input = _mm256_set1_ps(1.0);
  const aocommon::Avx256::VectorComplexFloat4 result{
      aocommon::Avx256::VectorFloat8{input}};

  BOOST_TEST(result[0] == (std::complex<float>{1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{1.0, 1.0}));
  BOOST_TEST(result[2] == (std::complex<float>{1.0, 1.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.0, 1.0}));
}

BOOST_AUTO_TEST_CASE(constructor_4_complex_float) {
  static_assert(std::is_nothrow_constructible_v<
                aocommon::Avx256::VectorComplexFloat4, std::complex<float>,
                std::complex<float>, std::complex<float>, std::complex<float>>);

  const aocommon::Avx256::VectorComplexFloat4 result{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
      std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_float_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat4,
                                      const std::complex<float>[4]>);
  const std::complex<float> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::VectorComplexFloat4 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_double_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::VectorComplexFloat4,
                                      const std::complex<double>[4]>);
  const std::complex<double> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::VectorComplexFloat4 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(operator_m256) {
  static_assert(
      noexcept(static_cast<__m256>(aocommon::Avx256::VectorComplexFloat4{
          static_cast<const std::complex<float>*>(nullptr)})));

  const aocommon::Avx256::VectorComplexFloat4 input{
      aocommon::Avx256::VectorFloat8{-1.0, 1.0, 3.75, -3.75, 99.0, -99.0, 1.5,
                                     -1.5}};

  BOOST_TEST(static_cast<__m256>(input)[0] == -1.0);
  BOOST_TEST(static_cast<__m256>(input)[1] == 1.0);
  BOOST_TEST(static_cast<__m256>(input)[2] == 3.75);
  BOOST_TEST(static_cast<__m256>(input)[3] == -3.75);
  BOOST_TEST(static_cast<__m256>(input)[4] == 99.0);
  BOOST_TEST(static_cast<__m256>(input)[5] == -99.0);
  BOOST_TEST(static_cast<__m256>(input)[6] == 1.5);
  BOOST_TEST(static_cast<__m256>(input)[7] == -1.5);
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(aocommon::Avx256::VectorComplexFloat4{
      static_cast<const std::complex<float>*>(nullptr)}
                             .Conjugate()));

  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {-10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {-10, 11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {-10, 11}, {-100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {-10, -11}, {-100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {-10, 11}, {-100, 101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::VectorComplexFloat4{
          {-1.0, -2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}
           .Conjugate()),
      (aocommon::Avx256::VectorComplexFloat4{
          {-1.0, 2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}
                         .Conjugate()),
                    (aocommon::Avx256::VectorComplexFloat4{
                        {-1.0, 2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}));
}

BOOST_AUTO_TEST_CASE(assign_to) {
  static_assert(noexcept(aocommon::Avx256::VectorComplexFloat4{
      static_cast<const std::complex<float>*>(nullptr)}
                             .AssignTo(nullptr)));

  const aocommon::Avx256::VectorComplexFloat4 input{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
      std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}};

  std::vector<std::complex<float>> result(6);
  input.AssignTo(std::addressof(result[1]));

  BOOST_TEST(result[0] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[2] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[3] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[4] == (std::complex<float>{1.5, -1.5}));
  BOOST_TEST(result[5] == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  aocommon::Avx256::VectorComplexFloat4 r{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::VectorComplexFloat4 value{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r[0] == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<float>{50.0, 55.0}));
  BOOST_TEST(r[2] == (std::complex<float>{500., 505.0}));
  BOOST_TEST(r[3] == (std::complex<float>{5000., 5005.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus_equal) {
  aocommon::Avx256::VectorComplexFloat4 r{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::VectorComplexFloat4 value{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  r -= value;
  static_assert(noexcept(r -= value));

  BOOST_TEST(r[0] == (std::complex<float>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<float>{-30, -33}));
  BOOST_TEST(r[2] == (std::complex<float>{-300, -303}));
  BOOST_TEST(r[3] == (std::complex<float>{-3000, -3003}));
}

BOOST_AUTO_TEST_CASE(operator_plus) {
  const aocommon::Avx256::VectorComplexFloat4 lhs{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::VectorComplexFloat4 rhs{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const aocommon::Avx256::VectorComplexFloat4 r = lhs + rhs;
  static_assert(noexcept(lhs + rhs));

  BOOST_TEST(r[0] == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<float>{50.0, 55.0}));
  BOOST_TEST(r[2] == (std::complex<float>{500., 505.0}));
  BOOST_TEST(r[3] == (std::complex<float>{5000., 5005.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus) {
  const aocommon::Avx256::VectorComplexFloat4 lhs{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::VectorComplexFloat4 rhs{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const aocommon::Avx256::VectorComplexFloat4 r = lhs - rhs;
  static_assert(noexcept(lhs - rhs));

  BOOST_TEST(r[0] == (std::complex<float>{-3, -6}));
  BOOST_TEST(r[1] == (std::complex<float>{-30, -33}));
  BOOST_TEST(r[2] == (std::complex<float>{-300, -303}));
  BOOST_TEST(r[3] == (std::complex<float>{-3000, -3003}));
}

BOOST_AUTO_TEST_CASE(multiply) {
  const aocommon::Avx256::VectorComplexFloat4 lhs{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::VectorComplexFloat4 rhs{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const aocommon::Avx256::VectorComplexFloat4 r = lhs * rhs;
  static_assert(noexcept(lhs * rhs));

  BOOST_TEST(r[0] == (std::complex<float>{-12.0, 16.0}));
  BOOST_TEST(r[1] == (std::complex<float>{-84.0, 880.0}));
  BOOST_TEST(r[2] == (std::complex<float>{-804.0, 80800.0}));
  BOOST_TEST(r[3] == (std::complex<float>{-8004.0, 8008000.0}));
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(aocommon::Avx256::VectorComplexFloat4{
                   static_cast<const std::complex<float>*>(nullptr)} ==
               aocommon::Avx256::VectorComplexFloat4{
                   static_cast<const std::complex<float>*>(nullptr)}));

  BOOST_TEST((aocommon::Avx256::VectorComplexFloat4{
                  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
              aocommon::Avx256::VectorComplexFloat4{
                  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}} ==
               aocommon::Avx256::VectorComplexFloat4{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_CASE(output) {
  const aocommon::Avx256::VectorComplexFloat4 input{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
      std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}};

  std::stringstream result;
  result << input;

  BOOST_CHECK_EQUAL(result.str(),
                    "[(-1,1), (3.75,-3.75), (99,-99), (1.5,-1.5)]");
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // defined(__AVX2__)
