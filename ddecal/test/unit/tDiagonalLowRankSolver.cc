// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/DiagonalLowRankSolver.h"

#include <vector>
#include <complex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor-blas/xlinalg.hpp>

BOOST_AUTO_TEST_SUITE(diagonal_low_rank_solver)

BOOST_AUTO_TEST_CASE(dominant_eigen_pair_2x2_real) {
  xt::xtensor<float, 2> m{{0.5f, 0.5f}, {0.2f, 0.8f}};
  xt::xtensor<std::complex<float>, 1> eigen_vector;
  const float eigen_value = dp3::ddecal::DominantEigenPair(m, eigen_vector, 20);
  BOOST_CHECK_CLOSE_FRACTION(eigen_value, 1.0f, 1e-6f);
  // The unscaled eigen vector is [1, 1]. Note that it may be multiplied by any
  // scalar, so we have to divide the ambiguity:
  eigen_vector = eigen_vector / eigen_vector[0];
  BOOST_CHECK_CLOSE_FRACTION(std::real(eigen_vector[1]), 1.0f, 1e-6f);
  BOOST_CHECK_LT(std::imag(eigen_vector[1]), 1e-6f);
}

BOOST_AUTO_TEST_CASE(dominant_eigen_pair_2x2_complex) {
  // Same test as above, but multiplied by (1 + 0.5i).
  xt::xtensor<std::complex<float>, 2> m{{0.5, 0.5}, {0.2, 0.8}};
  m = m * std::complex<float>(1.0, 0.5);
  xt::xtensor<std::complex<float>, 1> eigen_vector;
  const float eigen_value = dp3::ddecal::DominantEigenPair(m, eigen_vector, 20);
  BOOST_CHECK_CLOSE_FRACTION(eigen_value, 1.0f, 1e-6f);
  eigen_vector = eigen_vector / eigen_vector[0];
  BOOST_CHECK_CLOSE_FRACTION(std::real(eigen_vector[1]), 1.0f, 1e-6f);
  BOOST_CHECK_LT(std::imag(eigen_vector[1]), 1e-6f);
}

BOOST_AUTO_TEST_CASE(dominant_eigen_pair_3x3_real) {
  xt::xtensor<std::complex<float>, 1> t({std::complex<float>(1.0, 0.0),
                                         std::complex<float>(2.0, 0.0),
                                         std::complex<float>(3.0, 0.0)});
  xt::xtensor<std::complex<float>, 2> m = xt::linalg::outer(t, xt::conj(t));
  xt::xtensor<std::complex<float>, 1> eigen_vector;
  const float eigen_value = dp3::ddecal::DominantEigenPair(m, eigen_vector, 20);
  eigen_vector = eigen_vector * std::sqrt(eigen_value);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[0]), 1.0f, 1e-6f);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[1]), 2.0f, 1e-6f);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[2]), 3.0f, 1e-6f);
}

BOOST_AUTO_TEST_CASE(dominant_eigen_pair_3x3_complex) {
  xt::xtensor<std::complex<float>, 1> t({std::complex<float>(1.0, 0.0),
                                         std::complex<float>(4.0, 3.0),
                                         std::complex<float>(3.0, -0.4)});
  xt::xtensor<std::complex<float>, 2> m = xt::linalg::outer(t, xt::conj(t));
  xt::xtensor<std::complex<float>, 1> eigen_vector;
  const float eigen_value = dp3::ddecal::DominantEigenPair(m, eigen_vector, 20);
  eigen_vector = eigen_vector * std::sqrt(eigen_value);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[0]), std::fabs(t[0]),
                             1e-6f);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[1]), std::fabs(t[1]),
                             1e-6f);
  BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[2]), std::fabs(t[2]),
                             1e-6f);
}

BOOST_AUTO_TEST_CASE(dominant_eigen_pair_large_complex) {
  const size_t n = 100;
  xt::xtensor<std::complex<float>, 1> t;
  t.resize({n});
  std::uniform_real_distribution<float> distribution(0, n);
  std::mt19937 gen;
  for (std::complex<float>& v : t) {
    v = {distribution(gen), distribution(gen) * 0.5f};
  }
  xt::xtensor<std::complex<float>, 2> m = xt::linalg::outer(t, xt::conj(t));
  xt::xtensor<std::complex<float>, 1> eigen_vector;
  const float eigen_value = dp3::ddecal::DominantEigenPair(m, eigen_vector, 20);
  eigen_vector = eigen_vector * std::sqrt(eigen_value);
  for (size_t i = 0; i != n; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(std::fabs(eigen_vector[i]), std::fabs(t[i]),
                               1e-6f);
  }
  eigen_vector *= t[0] / eigen_vector[0];
  for (size_t i = 0; i != n; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(eigen_vector[i].real(), t[i].real(), 1e-5f);
  }
}

BOOST_AUTO_TEST_SUITE_END()
