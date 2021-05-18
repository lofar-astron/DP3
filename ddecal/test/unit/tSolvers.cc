// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTester.h"

#include "../../gain_solvers/DiagonalSolver.h"
#include "../../gain_solvers/FullJonesSolver.h"
#include "../../gain_solvers/HybridSolver.h"
#include "../../gain_solvers/IterativeDiagonalSolver.h"
#include "../../gain_solvers/IterativeScalarSolver.h"
#include "../../gain_solvers/ScalarSolver.h"
#include "../../gain_solvers/SolverBuffer.h"

#include "../../linear_solvers/LSMRSolver.h"
#include "../../linear_solvers/NormalEquationsSolver.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::LLSSolverType;
using dp3::ddecal::SolverBuffer;
using dp3::ddecal::test::SolverTester;

// The solvers test suite also contains tests that run using a separate ctest
// test, since they take much time. These tests have the 'slow' label.
BOOST_AUTO_TEST_SUITE(solvers)

#ifdef USE_LSMR
BOOST_AUTO_TEST_CASE(lsmr_solver) {
  int n = 2;
  int m = 4;
  std::vector<std::complex<float>> A{1.0, 1.5, 3.5, 2.0,
                                     0.5, 0.7, 0.8, 0.4};  // column major
  std::vector<std::complex<float>> b{1.0, 2.0, 1.5, 1.2};

  dp3::ddecal::LSMRSolver solver(m, n, 1);
  solver.Solve(A.data(), b.data());

  BOOST_CHECK_CLOSE(b[0].real(), -0.14141126, 1.0E-3);
  BOOST_CHECK_CLOSE(b[1].real(), 2.79757662, 1.0E-3);
}

BOOST_FIXTURE_TEST_CASE(scalar_solver_lsmr, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::ScalarSolver solver;
  InitializeSolver(solver);
  solver.SetAccuracy(1e-9);
  solver.SetLLSSolverType(LLSSolverType::LSMR, 1.0E-2, 1.0E-2);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}
#endif

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

BOOST_FIXTURE_TEST_CASE(scalar_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::ScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR, 0.0, 0.0);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_scalar_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeScalarSolver solver;
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0e-3);
}

BOOST_FIXTURE_TEST_CASE(scalar_solver_normaleq, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::ScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::NORMAL_EQUATIONS, 0.0, 0.0);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

#ifdef USE_LSMR
BOOST_FIXTURE_TEST_CASE(diagonal_solver_lsmr, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::DiagonalSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::LSMR, 1.0E-7, 1.0E-2);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetDiagonalSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(2e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}
#endif

BOOST_FIXTURE_TEST_CASE(diagonal_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::DiagonalSolver solver;
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetDiagonalSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(2e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_diagonal_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeDiagonalSolver solver;
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetDiagonalSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(hybrid_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  auto direction_solver = boost::make_unique<dp3::ddecal::DiagonalSolver>();
  InitializeSolver(*direction_solver);
  direction_solver->SetMaxIterations(kMaxIterations / 10);

  auto iterative_solver =
      boost::make_unique<dp3::ddecal::IterativeDiagonalSolver>();
  InitializeSolver(*iterative_solver);

  dp3::ddecal::HybridSolver solver;
  solver.AddSolver(std::move(iterative_solver));
  solver.AddSolver(std::move(direction_solver));
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 2u);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[0], &solver);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[1], &solver);

  SetDiagonalSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(full_jones_solver, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::FullJonesSolver solver;
  InitializeSolver(solver);
  solver.AddConstraint(boost::make_unique<dp3::ddecal::DiagonalConstraint>(4));

  SetDiagonalSolutions();

  const SolverBuffer& solver_buffer = FillData();

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 4u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  dp3::ddecal::SolverBase::SolveResult result;

  // The full jones test uses full matrices as solutions and copies the
  // diagonals into the solver solutions from the SolverTester fixture. This
  // way, the test can reuse SetDiagonalSolutions() and CheckDiagonalResults().
  std::vector<std::vector<std::complex<double>>> solutions(kNChannelBlocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(kNDirections * kNAntennas * 4, 0.0);
    for (size_t i = 0; i != kNDirections * kNAntennas * 4; i += 4) {
      vec[i] = 1.0;
      vec[i + 3] = 1.0;
    }
  }

  result = solver.Solve(solver_buffer, solutions, 0.0, nullptr);

  // Convert full matrices to diagonals
  std::vector<std::vector<std::complex<double>>>& diagonals =
      GetSolverSolutions();
  for (size_t chBlock = 0; chBlock != solutions.size(); ++chBlock) {
    for (size_t s = 0; s != solutions[chBlock].size() / 4; s++) {
      diagonals[chBlock][s * 2] = solutions[chBlock][s * 4];
      diagonals[chBlock][s * 2 + 1] = solutions[chBlock][s * 4 + 3];
    }
  }

  CheckDiagonalResults(2e-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(min_iterations, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::IterativeScalarSolver solver;
  InitializeSolver(solver);
  solver.SetMinIterations(10);
  // very large tolerance on purpose to stop directly once min iters are reached
  solver.SetAccuracy(1e8);

  SetScalarSolutions();

  const SolverBuffer& solver_buffer = FillData();

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(solver_buffer, GetSolverSolutions(), 0.0, nullptr);
  BOOST_CHECK_EQUAL(result.iterations, 10U);
}

BOOST_AUTO_TEST_SUITE_END()
