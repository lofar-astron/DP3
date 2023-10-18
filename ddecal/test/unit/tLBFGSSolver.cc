// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/LBFGSSolver.h"

#include <boost/test/unit_test.hpp>

#ifdef HAVE_LIBDIRAC
using dp3::ddecal::LBFGSSolver;
using dp3::ddecal::SolutionTensor;

BOOST_AUTO_TEST_SUITE(lbfgs_solver)

BOOST_AUTO_TEST_CASE(split_solutions) {
  {
    const std::vector<std::complex<double>> empty_input;
    const xt::xtensor<double, 1> output =
        LBFGSSolver::SplitSolutions(empty_input);
    BOOST_TEST(output.size() == 0u);
  }
  {
    const double kReal = 42.0;
    const double kImaginary = 43.0;
    const std::vector<std::complex<double>> single_input{{kReal, kImaginary}};

    const xt::xtensor<double, 1> output =
        LBFGSSolver::SplitSolutions(single_input);

    BOOST_TEST(output.size() == 2u);
    BOOST_TEST(output[0] == kReal);
    BOOST_TEST(output[1] == kImaginary);
  }
  {
    const double kFirstReal = 42.0;
    const double kFirstImaginary = 142.0;
    const std::vector<std::complex<double>> multiple_inputs{
        {kFirstReal + 0.0, kFirstImaginary + 0.0},
        {kFirstReal + 1.0, kFirstImaginary + 1.0},
        {kFirstReal + 2.0, kFirstImaginary + 2.0},
        {kFirstReal + 3.0, kFirstImaginary + 3.0},
        {kFirstReal + 4.0, kFirstImaginary + 4.0}};

    const xt::xtensor<double, 1> output =
        LBFGSSolver::SplitSolutions(multiple_inputs);

    BOOST_TEST_REQUIRE(output.size() == 10u);
    for (int i = 0; i < 5; ++i) {
      BOOST_TEST(output[i] == kFirstReal + i);
      BOOST_TEST(output[5 + i] == kFirstImaginary + i);
    }
  }
}

BOOST_AUTO_TEST_CASE(merge_solutions) {
  {
    const xt::xtensor<double, 1> empty_input;
    SolutionTensor empty_solution;

    BOOST_CHECK_NO_THROW(
        LBFGSSolver::MergeSolutions(empty_solution, 42, empty_input));
  }
  {  // Test with a single complex number. Use multiple channel blocks.
    const double kReal = 42.0;
    const double kImaginary = 43.0;
    const std::complex<double> kDefaultValue{0.0, 0.0};
    const std::complex<double> kNewValue{kReal, kImaginary};
    const int kChannelBlock = 41;
    const size_t kNChannelBlocks = 42;
    const xt::xtensor<double, 1> d_storage{kReal, kImaginary};
    const std::array<size_t, 4> shape{kNChannelBlocks, 1, 1, 1};
    SolutionTensor solution(shape, kDefaultValue);

    LBFGSSolver::MergeSolutions(solution, kChannelBlock, d_storage);

    for (int i = 0; i < kChannelBlock; ++i) {
      BOOST_TEST(solution(i, 0, 0, 0) == kDefaultValue);
    }
    BOOST_TEST(solution(kChannelBlock, 0, 0, 0) == kNewValue);
  }
  {  // Test with five complex numbers. Use a single channel block.
    const double kFirstReal = 42.0;
    const double kFirstImaginary = 142.0;
    const size_t kChannelBlock = 0;
    const int kNValues = 5;
    const std::array<size_t, 1> d_storage_shape{kNValues * 2};
    xt::xtensor<double, 1> d_storage(d_storage_shape);
    for (int i = 0; i < kNValues; ++i) {
      d_storage(i) = kFirstReal + i;
      d_storage(i + kNValues) = kFirstImaginary + i;
    }
    const std::array<size_t, 4> shape{1, kNValues, 1, 1};
    SolutionTensor solution(shape);

    LBFGSSolver::MergeSolutions(solution, kChannelBlock, d_storage);

    for (int i = 0; i < kNValues; ++i) {
      const std::complex<double> expected_value{kFirstReal + i,
                                                kFirstImaginary + i};
      BOOST_TEST(solution(0, i, 0, 0) == expected_value);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
#endif  // HAVE_LIBDIRAC
