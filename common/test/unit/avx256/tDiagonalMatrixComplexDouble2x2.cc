// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/DiagonalMatrixComplexDouble2x2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

#if defined(__AVX2__)

// Silences the diagnostic "ignoring attributes on template argument ‘__m256d’"
// This diagnostic is issued on code that's part of libstdc++.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

// operator[] is tested in other tests.

BOOST_AUTO_TEST_SUITE(DiagonalMatrixComplexDouble2x2)

static_assert(!std::is_default_constructible_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);
static_assert(std::is_nothrow_destructible_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);
static_assert(std::is_nothrow_copy_constructible_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);
static_assert(std::is_nothrow_move_constructible_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);
static_assert(std::is_nothrow_copy_assignable_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);
static_assert(std::is_nothrow_move_assignable_v<
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2>);

BOOST_AUTO_TEST_CASE(constructor_2_complex_double) {
  static_assert(std::is_nothrow_constructible_v<
                aocommon::Avx256::DiagonalMatrixComplexDouble2x2,
                std::complex<double>, std::complex<double>>);

  const aocommon::Avx256::DiagonalMatrixComplexDouble2x2 result{
      std::complex<double>{-1.0, 1.0}, std::complex<double>{3.75, -3.75}};

  BOOST_TEST(result[0] == (std::complex<double>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<double>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_double_pointer) {
  static_assert(std::is_nothrow_constructible_v<
                aocommon::Avx256::DiagonalMatrixComplexDouble2x2,
                const std::complex<double>[4]>);
  const std::complex<double> input[] = {std::complex<double>{-1.0, 1.0},
                                        std::complex<double>{3.75, -3.75}};
  const aocommon::Avx256::DiagonalMatrixComplexDouble2x2 result{input};

  BOOST_TEST((result[0] == std::complex<double>{-1.0, 1.0}));
  BOOST_TEST((result[1] == std::complex<double>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
      static_cast<const std::complex<double>*>(nullptr)}
                             .Conjugate()));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{1.0, 2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{1.0, -2.0},
                                                        {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0},
                                                        {10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .Conjugate()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .Conjugate()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .Conjugate()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .Conjugate()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(hermitian_transpose) {
  static_assert(noexcept(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
      static_cast<const std::complex<double>*>(nullptr)}
                             .HermitianTranspose()));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{1.0, 2.0}, {10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{1.0, -2.0},
                                                        {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0}, {10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0},
                                                        {10, -11}}));
  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {10, -11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .HermitianTranspose()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .HermitianTranspose()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL(
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, -2.0}, {-10, 11}}
           .HermitianTranspose()),
      (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{-1.0, 2.0},
                                                        {-10, -11}}));
  BOOST_CHECK_EQUAL((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, -2.0}, {-10, -11}}
                         .HermitianTranspose()),
                    (aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                        {-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                   static_cast<const std::complex<double>*>(nullptr)} ==
               aocommon::Avx256::DiagonalMatrixComplexDouble2x2{
                   static_cast<const std::complex<double>*>(nullptr)}));

  BOOST_TEST((aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                               {0.0, 0.0}} ==
              aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                               {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{42.0, 0.0},
                                                                {0.0, 0.0}} ==
               aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 42.0},
                                                                {0.0, 0.0}} ==
               aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {42.0, 0.0}} ==
               aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {0.0, 0.0}}));
  BOOST_TEST(!(aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {0.0, 42.0}} ==
               aocommon::Avx256::DiagonalMatrixComplexDouble2x2{{0.0, 0.0},
                                                                {0.0, 0.0}}));
}

BOOST_AUTO_TEST_SUITE_END()

#pragma GCC diagnostic pop

#endif  // defined(__AVX2__)
