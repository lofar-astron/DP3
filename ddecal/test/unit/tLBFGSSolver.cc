// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/LBFGSSolver.h"

#include <boost/test/unit_test.hpp>

#include "SolverTester.h"

using dp3::ddecal::LBFGSSolver;
using dp3::ddecal::SolutionTensor;
using dp3::ddecal::SolveData;
using dp3::ddecal::test::SolverTester;

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

BOOST_FIXTURE_TEST_CASE(diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches, 0.0, 0.0,
                                  dp3::ddecal::LBFGSSolver::kDiagonal);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);
  }

  CheckDiagonalResults(2.0E-2);
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches, 0.0, 0.0,
                                  dp3::ddecal::LBFGSSolver::kScalar);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);
  }

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(full_jones, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  // Since we have more unknowns, reduce the minibatch size
  // to increase the constraints per-minibatch
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations * 2, kHistory,
                                  kMinibatches / 2, 0.0, 0.0,
                                  dp3::ddecal::LBFGSSolver::kFull);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  // The full jones test uses full matrices as solutions and copies the
  // diagonals into the solver solutions from the SolverTester fixture. This
  // way, the test can reuse SetDiagonalSolutions() and CheckDiagonalResults().
  std::vector<std::vector<std::complex<double>>> solutions(kNChannelBlocks);

  // Initialize unit-matrices as initial values
  for (auto& solution : solutions) {
    solution.assign(NSolutions() * kNAntennas * 4, 0.0);
    for (size_t i = 0; i != NSolutions() * kNAntennas * 4; i += 4) {
      solution[i] = 1.0;
      solution[i + 3] = 1.0;
    }
  }

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, solutions, 0.0, nullptr);
  }

  // Convert full matrices to diagonals
  std::vector<std::vector<std::complex<double>>>& diagonals =
      GetSolverSolutions();
  for (size_t ch_block = 0; ch_block != solutions.size(); ++ch_block) {
    for (size_t s = 0; s != solutions[ch_block].size() / 4; ++s) {
      diagonals[ch_block][s * 2] = solutions[ch_block][s * 4];
      diagonals[ch_block][s * 2 + 1] = solutions[ch_block][s * 4 + 3];
    }
  }

  CheckDiagonalResults(2.0e-2);
  // The solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(bounded_diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations * 2, kHistory,
                                  kMinibatches / 2, kMinSolution, kMaxSolution,
                                  dp3::ddecal::LBFGSSolver::kDiagonal);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);
  }

  CheckDiagonalResults(10.0);
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(bounded_scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches, kMinSolution, kMaxSolution,
                                  dp3::ddecal::LBFGSSolver::kScalar);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);
  }

  CheckScalarResults(10.0);
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(bounded_full_jones, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  // Since we have more unknowns, reduce the minibatch size
  // to increase the constraints per-minibatch
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations * 2, kHistory,
                                  kMinibatches / 2, kMinSolution, kMaxSolution,
                                  dp3::ddecal::LBFGSSolver::kFull);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNAntennas,
                       std::vector<size_t>(kNDirections, 1), Antennas1(),
                       Antennas2(), false);

  // The full jones test uses full matrices as solutions and copies the
  // diagonals into the solver solutions from the SolverTester fixture. This
  // way, the test can reuse SetDiagonalSolutions() and CheckDiagonalResults().
  std::vector<std::vector<std::complex<double>>> solutions(kNChannelBlocks);

  // Initialize unit-matrices as initial values
  for (auto& solution : solutions) {
    solution.assign(NSolutions() * kNAntennas * 4, 0.0);
    for (size_t i = 0; i != NSolutions() * kNAntennas * 4; i += 4) {
      solution[i] = 1.0;
      solution[i + 3] = 1.0;
    }
  }

  dp3::ddecal::SolverBase::SolveResult result;
  for (size_t n_epoch = 0; n_epoch < kEpochs; n_epoch++) {
    result = solver.Solve(data, solutions, 0.0, nullptr);
  }

  // Convert full matrices to diagonals
  std::vector<std::vector<std::complex<double>>>& diagonals =
      GetSolverSolutions();
  for (size_t ch_block = 0; ch_block != solutions.size(); ++ch_block) {
    for (size_t s = 0; s != solutions[ch_block].size() / 4; ++s) {
      diagonals[ch_block][s * 2] = solutions[ch_block][s * 4];
      diagonals[ch_block][s * 2 + 1] = solutions[ch_block][s * 4 + 3];
    }
  }

  CheckDiagonalResults(10.0);
  // The solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_AUTO_TEST_SUITE_END()
