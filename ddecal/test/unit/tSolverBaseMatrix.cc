// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ddecal/gain_solvers/SolverBase.h"

#include <boost/test/unit_test.hpp>

// There are no tests for the default constuctor since the operations are
// either undefined behaviour or have unspecified results.
// The same holds true for testing with 0 rows or 0 columns.

BOOST_AUTO_TEST_SUITE(solver_base_matrix, *boost::unit_test::tolerance(0.0f))

BOOST_AUTO_TEST_CASE(constructor) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_CASE(set_zero) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(1.0, 1.0);

  matrix.SetZero();

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_CASE(operator_parenthesis) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(column, row);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(column, row));
}

BOOST_AUTO_TEST_CASE(data) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(column, row);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      BOOST_CHECK(matrix.data()[column + 3 * row] ==
                  dp3::ddecal::SolverBase::Complex(column, row));
}

BOOST_AUTO_TEST_CASE(reset_shrink) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(1.0, 1.0);

  matrix.Reset(1, 1);

  BOOST_CHECK(matrix(0, 0) == dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_CASE(reset_equal_column_row) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(1.0, 1.0);

  matrix.Reset(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_CASE(reset_swap_row_column) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(1.0, 1.0);

  matrix.Reset(5, 3);

  for (int column = 0; column < 5; ++column)
    for (int row = 0; row < 3; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_CASE(reset_grow) {
  dp3::ddecal::SolverBase::Matrix matrix(3, 5);

  for (int column = 0; column < 3; ++column)
    for (int row = 0; row < 5; ++row)
      matrix(column, row) = dp3::ddecal::SolverBase::Complex(1.0, 1.0);

  matrix.Reset(10, 10);

  for (int column = 0; column < 10; ++column)
    for (int row = 0; row < 10; ++row)
      BOOST_CHECK(matrix(column, row) ==
                  dp3::ddecal::SolverBase::Complex(0.0, 0.0));
}

BOOST_AUTO_TEST_SUITE_END()
