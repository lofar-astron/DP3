// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTester.h"

#include "../../gain_solvers/DiagonalLowRankSolver.h"
#include "../../gain_solvers/DiagonalSolver.h"
#include "../../gain_solvers/FullJonesSolver.h"
#include "../../gain_solvers/HybridSolver.h"
#include "../../gain_solvers/IterativeDiagonalSolver.h"
#if defined(HAVE_CUDA_SOLVER)
#include "../../gain_solvers/IterativeDiagonalSolverCuda.h"
#endif
#include "../../gain_solvers/IterativeFullJonesSolver.h"
#include "../../gain_solvers/IterativeScalarSolver.h"
#include "../../gain_solvers/ScalarSolver.h"
#ifdef HAVE_LIBDIRAC
#include "../../gain_solvers/LBFGSSolver.h"
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::LLSSolverType;
using dp3::ddecal::SolveData;
using dp3::ddecal::test::SolverTester;

// The solver test suite contains tests that run using a separate
// ctest test, since they take much time. These tests have the 'slow' label.
BOOST_AUTO_TEST_SUITE(solvers)

BOOST_FIXTURE_TEST_CASE(diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::DiagonalSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(2.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::ScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(diagonal_low_rank_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  constexpr size_t kSolveIterations = 100;
  constexpr size_t kNLraIterations = 25;
  constexpr size_t kNPowerIterations = 10;
  dp3::ddecal::DiagonalLowRankSolver solver(kNLraIterations, kNPowerIterations);
  InitializeSolver(solver);
  solver.SetMaxIterations(kSolveIterations);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), true);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1.0E-3);
  BOOST_CHECK_EQUAL(result.iterations, kSolveIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(diagonal_low_rank_full_step_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  constexpr size_t kSolveIterations = 25;
  constexpr size_t kNLraIterations = 25;
  constexpr size_t kNPowerIterations = 10;
  dp3::ddecal::DiagonalLowRankSolver solver(kNLraIterations, kNPowerIterations);
  InitializeSolver(solver);
  solver.SetStepSize(1.0);
  solver.SetMaxIterations(kSolveIterations);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), true);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1.0E-3);
  BOOST_CHECK_EQUAL(result.iterations, kSolveIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::IterativeScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  // The iterative solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_scalar_dd_intervals, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(true);
  dp3::ddecal::IterativeScalarSolver solver;
  InitializeSolver(solver);

  const std::vector<dp3::base::DPBuffer> data_buffers = FillDdIntervalData();
  const SolveData data(data_buffers, CreateDirectionNames(), kNChannelBlocks,
                       kNAntennas, NSolutionsPerDirection(), Antennas1(),
                       Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0e-3);
}

BOOST_FIXTURE_TEST_CASE(hybrid, SolverTester,
                        *boost::unit_test::label("slow")) {
  auto direction_solver = std::make_unique<dp3::ddecal::ScalarSolver>();
  direction_solver->SetMaxIterations(kMaxIterations / 10);

  auto iterative_solver =
      std::make_unique<dp3::ddecal::IterativeScalarSolver>();
  iterative_solver->SetMaxIterations(kMaxIterations);

  SetScalarSolutions(false);
  dp3::ddecal::HybridSolver solver;
  solver.AddSolver(std::move(iterative_solver));
  solver.AddSolver(std::move(direction_solver));
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 2u);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[0], &solver);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[1], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

inline void TestIterativeDiagonal(SolverTester& solver_tester,
                                  dp3::ddecal::SolverBase& solver) {
  solver_tester.SetDiagonalSolutions(false);
  solver_tester.InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer =
      solver_tester.FillBDAData();
  const SolveData data(solver_buffer, SolverTester::kNChannelBlocks,
                       SolverTester::kNDirections, SolverTester::kNAntennas,
                       solver_tester.Antennas1(), solver_tester.Antennas2(),
                       false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, solver_tester.GetSolverSolutions(), 0.0, nullptr);

  solver_tester.CheckDiagonalResults(1.0E-2);
  // The iterative solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, SolverTester::kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeDiagonalSolver solver;
  TestIterativeDiagonal(*this, solver);
}

#if defined(HAVE_CUDA_SOLVER)
BOOST_FIXTURE_TEST_CASE(iterative_diagonal_cuda, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeDiagonalSolverCuda solver;
  TestIterativeDiagonal(*this, solver);
}
BOOST_FIXTURE_TEST_CASE(iterative_diagonal_cuda_keep_buffers, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeDiagonalSolverCuda solver{true};
  TestIterativeDiagonal(*this, solver);
}

#endif

BOOST_FIXTURE_TEST_CASE(iterative_diagonal_dd_intervals, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(true);
  dp3::ddecal::IterativeDiagonalSolver solver;
  InitializeSolver(solver);

  const std::vector<dp3::base::DPBuffer> data_buffers = FillDdIntervalData();
  const SolveData data(data_buffers, CreateDirectionNames(), kNChannelBlocks,
                       kNAntennas, NSolutionsPerDirection(), Antennas1(),
                       Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1.0e-2);
}

BOOST_FIXTURE_TEST_CASE(full_jones, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::FullJonesSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

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

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, solutions, 0.0, nullptr);

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

BOOST_FIXTURE_TEST_CASE(iterative_full_jones, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::IterativeFullJonesSolver solver;
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2(), false);

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

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, solutions, 0.0, nullptr);

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

BOOST_FIXTURE_TEST_CASE(iterative_full_jones_dd_intervals, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(true);
  dp3::ddecal::IterativeFullJonesSolver solver;
  InitializeSolver(solver);

  const std::vector<dp3::base::DPBuffer> data_buffers = FillDdIntervalData();
  const SolveData data(data_buffers, CreateDirectionNames(), kNChannelBlocks,
                       kNAntennas, NSolutionsPerDirection(), Antennas1(),
                       Antennas2());

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

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, solutions, 0.0, nullptr);

  // Convert full matrices to diagonals
  std::vector<std::vector<std::complex<double>>>& diagonals =
      GetSolverSolutions();
  for (size_t ch_block = 0; ch_block != solutions.size(); ++ch_block) {
    for (size_t s = 0; s != solutions[ch_block].size() / 4; ++s) {
      diagonals[ch_block][s * 2] = solutions[ch_block][s * 4];
      diagonals[ch_block][s * 2 + 1] = solutions[ch_block][s * 4 + 3];
    }
  }

  CheckDiagonalResults(1.0e-2);
}

BOOST_FIXTURE_TEST_CASE(scalar_normaleq, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::ScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::NORMAL_EQUATIONS);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(min_iterations, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::IterativeScalarSolver solver;
  InitializeSolver(solver);
  solver.SetMinIterations(10);
  // very large tolerance on purpose to stop directly once min iters are reached
  solver.SetAccuracy(1e8);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);
  BOOST_CHECK_EQUAL(result.iterations, 10U);
}

#ifdef HAVE_LIBDIRAC
BOOST_FIXTURE_TEST_CASE(lbfgs_diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches,
                                  dp3::ddecal::LBFGSSolver::kDiagonal);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(2.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(lbfgs_scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetScalarSolutions(false);
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches,
                                  dp3::ddecal::LBFGSSolver::kScalar);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(lbfgs_full_jones, SolverTester,
                        *boost::unit_test::label("slow")) {
  SetDiagonalSolutions(false);
  // Since we have more unknowns, reduce the minibatch size
  // to increase the constraints per-minibatch
  dp3::ddecal::LBFGSSolver solver(kRobustDOF, kBatchIterations, kHistory,
                                  kMinibatches / 2,
                                  dp3::ddecal::LBFGSSolver::kFull);
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  const SolveData data(solver_buffer, kNChannelBlocks, kNDirections, kNAntennas,
                       Antennas1(), Antennas2(), false);

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

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, solutions, 0.0, nullptr);

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
#endif /* HAVE_LIBDIRAC */

BOOST_AUTO_TEST_SUITE_END()
