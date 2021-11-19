// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../linear_solvers/NormalEquationsSolver.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

BOOST_AUTO_TEST_SUITE(linear_solvers)

BOOST_AUTO_TEST_CASE(normaleq_solver) {
  int n = 2;
  int m = 4;
  std::vector<std::complex<float>> A{1.0, 1.5, 3.5, 2.0,
                                     0.5, 0.7, 0.8, 0.4};  // column major
  std::vector<std::complex<float>> b{1.0, 2.0, 1.5, 1.2};

  dp3::ddecal::NormalEquationsSolver solver(m, n, 1);
  solver.Solve(A.data(), b.data());

  BOOST_CHECK_CLOSE(b[0].real(), -0.14141126, 1.0E-3);
  BOOST_CHECK_CLOSE(b[1].real(), 2.79757662, 1.0E-3);
}

BOOST_AUTO_TEST_SUITE_END()
