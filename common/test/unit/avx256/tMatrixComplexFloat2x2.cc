// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/MatrixComplexFloat2x2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(MaxtrixComplexFloat2x2)

static_assert(
    !std::is_default_constructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2>);
static_assert(
    std::is_nothrow_destructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2>);
static_assert(std::is_nothrow_copy_constructible_v<
              aocommon::Avx256::MaxtrixComplexFloat2x2>);
static_assert(std::is_nothrow_move_constructible_v<
              aocommon::Avx256::MaxtrixComplexFloat2x2>);
static_assert(std::is_nothrow_copy_assignable_v<
              aocommon::Avx256::MaxtrixComplexFloat2x2>);
static_assert(std::is_nothrow_move_assignable_v<
              aocommon::Avx256::MaxtrixComplexFloat2x2>);

BOOST_AUTO_TEST_CASE(constructor_vector_complex_float_4) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2,
                                      aocommon::Avx256::VectorComplexFloat4>);

  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{
      aocommon::Avx256::VectorComplexFloat4{
          std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
          std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}}};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_4_complex_floats) {
  static_assert(std::is_nothrow_constructible_v<
                aocommon::Avx256::MaxtrixComplexFloat2x2, std::complex<float>,
                std::complex<float>, std::complex<float>, std::complex<float>>);

  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
      std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_float_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2,
                                      const std::complex<float>[4]>);
  const std::complex<float> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_double_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2,
                                      const std::complex<double>[4]>);
  const std::complex<double> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_MC2x2F) {
  static_assert(
      std::is_nothrow_constructible_v<aocommon::Avx256::MaxtrixComplexFloat2x2,
                                      const aocommon::MC2x2F&>);
  const aocommon::MC2x2F input{
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{input};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[2] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[3] == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(operator_MC2x2F) {
  const aocommon::MC2x2F input{
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const aocommon::Avx256::MaxtrixComplexFloat2x2 result{input};

  static_assert(noexcept(static_cast<aocommon::MC2x2F>(result)));

  BOOST_TEST(input[0] == static_cast<aocommon::MC2x2F>(result)[0]);
  BOOST_TEST(input[1] == static_cast<aocommon::MC2x2F>(result)[1]);
  BOOST_TEST(input[2] == static_cast<aocommon::MC2x2F>(result)[2]);
  BOOST_TEST(input[3] == static_cast<aocommon::MC2x2F>(result)[3]);
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(aocommon::Avx256::MaxtrixComplexFloat2x2{
      static_cast<const std::complex<float>*>(nullptr)}
                             .Conjugate()));

  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {-10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {-10, 11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {-10, 11}, {-100, 101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {-10, -11}, {-100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {-10, 11}, {-100, 101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::MaxtrixComplexFloat2x2{
          {-1.0, -2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}
           .Conjugate()),
      (aocommon::Avx256::MaxtrixComplexFloat2x2{
          {-1.0, 2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}
                         .Conjugate()),
                    (aocommon::Avx256::MaxtrixComplexFloat2x2{
                        {-1.0, 2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}));
}

BOOST_AUTO_TEST_CASE(transpose) {
  static_assert(noexcept(aocommon::Avx256::MaxtrixComplexFloat2x2{
      static_cast<const std::complex<float>*>(nullptr)}
                             .Transpose()));
  const aocommon::MC2x2F input{
      std::complex<float>{1.0, 2.0}, std::complex<float>{10, 11},
      std::complex<float>{100, 101}, std::complex<float>{1000, 1001}};

  const aocommon::Avx256::MaxtrixComplexFloat2x2 expected{input.Transpose()};

  const aocommon::Avx256::MaxtrixComplexFloat2x2 result =
      aocommon::Avx256::MaxtrixComplexFloat2x2{input}.Transpose();

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(hermitian_transpose) {
  static_assert(noexcept(aocommon::Avx256::MaxtrixComplexFloat2x2{
      static_cast<const std::complex<float>*>(nullptr)}
                             .HermitianTranspose()));
  const std::vector<aocommon::MC2x2F> inputs{
      aocommon::MC2x2F{
          std::complex<float>{1.0, 2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, 2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, 101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{-1000, 1001}},
      aocommon::MC2x2F{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{-1000, -1001}}};

  for (const auto& input : inputs) {
    const aocommon::Avx256::MaxtrixComplexFloat2x2 expected{
        input.HermTranspose()};

    const aocommon::Avx256::MaxtrixComplexFloat2x2 result =
        aocommon::Avx256::MaxtrixComplexFloat2x2{input}.HermitianTranspose();

    BOOST_CHECK_EQUAL(result, expected);
  }
}

BOOST_AUTO_TEST_CASE(multiply) {
  const aocommon::Avx256::MaxtrixComplexFloat2x2 lhs{
      {1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::Avx256::MaxtrixComplexFloat2x2 rhs{
      {4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  static_assert(noexcept(lhs * rhs));
  const aocommon::Avx256::MaxtrixComplexFloat2x2 r = lhs * rhs;

  BOOST_CHECK_CLOSE(r[0].real(), -456, 1e-6);
  BOOST_CHECK_CLOSE(r[0].imag(), 8456, 1e-6);
  BOOST_CHECK_CLOSE(r[1].real(), -4092, 1e-6);
  BOOST_CHECK_CLOSE(r[1].imag(), 84164, 1e-6);
  BOOST_CHECK_CLOSE(r[2].real(), -4812, 1e-6);
  BOOST_CHECK_CLOSE(r[2].imag(), 805604, 1e-6);
  BOOST_CHECK_CLOSE(r[3].real(), -8448, 1e-6);
  BOOST_CHECK_CLOSE(r[3].imag(), 8016440, 1e-6);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   static_cast<const std::complex<float>*>(nullptr)} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   static_cast<const std::complex<float>*>(nullptr)}));

  BOOST_TEST((aocommon::Avx256::MaxtrixComplexFloat2x2{
                  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
              aocommon::Avx256::MaxtrixComplexFloat2x2{
                  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}} ==
               aocommon::Avx256::MaxtrixComplexFloat2x2{
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_CASE(output) {
  const aocommon::Avx256::MaxtrixComplexFloat2x2 input{
      std::complex<float>{-1.0, 1.0}, std::complex<float>{3.75, -3.75},
      std::complex<float>{99.0, -99.0}, std::complex<float>{1.5, -1.5}};

  std::stringstream result;
  result << input;

  BOOST_CHECK_EQUAL(result.str(),
                    "[{(-1,1), (3.75,-3.75)}, {(99,-99), (1.5,-1.5)}]");
}

BOOST_AUTO_TEST_CASE(herm_transpose) {
  static_assert(noexcept(HermTranspose(aocommon::Avx256::MaxtrixComplexFloat2x2{
      static_cast<const std::complex<float>*>(nullptr)})));
  const std::vector<aocommon::Avx256::MaxtrixComplexFloat2x2> inputs{
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{1.0, 2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, 2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, 11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, 101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{-1000, 1001}},
      aocommon::Avx256::MaxtrixComplexFloat2x2{
          std::complex<float>{-1.0, -2.0}, std::complex<float>{-10, -11},
          std::complex<float>{-100, -101}, std::complex<float>{-1000, -1001}}};

  for (const auto& input : inputs)
    BOOST_CHECK_EQUAL(input.HermitianTranspose(), HermTranspose(input));
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // defined(__AVX2__)
