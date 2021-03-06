// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTester.h"

#include "../../gain_solvers/BdaDiagonalSolver.h"
#include "../../gain_solvers/BdaHybridSolver.h"
#include "../../gain_solvers/BdaIterativeDiagonalSolver.h"
#include "../../gain_solvers/BdaIterativeScalarSolver.h"
#include "../../gain_solvers/BdaScalarSolver.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/make_unique.hpp>

using dp3::ddecal::LLSSolverType;
using dp3::ddecal::test::SolverTester;

// The BDA solvers test suite also contains tests that run using a separate
// ctest test, since they take much time. These tests have the 'slow' label.
BOOST_AUTO_TEST_SUITE(bda_solvers)

BOOST_FIXTURE_TEST_CASE(diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::BdaDiagonalSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR, 0.0, 0.0);

  SetDiagonalSolutions();

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(2.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::BdaScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR, 0.0, 0.0);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_scalar, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::BdaIterativeScalarSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR, 0.0, 0.0);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetScalarSolutions();

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  // The iterative solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(hybrid, SolverTester,
                        *boost::unit_test::label("slow")) {
  auto direction_solver = boost::make_unique<dp3::ddecal::BdaScalarSolver>();
  direction_solver->SetMaxIterations(kMaxIterations / 10);

  auto iterative_solver =
      boost::make_unique<dp3::ddecal::BdaIterativeScalarSolver>();
  iterative_solver->SetMaxIterations(kMaxIterations);

  dp3::ddecal::BdaHybridSolver solver;
  solver.AddSolver(std::move(iterative_solver));
  solver.AddSolver(std::move(direction_solver));
  InitializeSolver(solver);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 1u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 2u);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[0], &solver);
  BOOST_CHECK_NE(solver.ConstraintSolvers()[1], &solver);

  SetScalarSolutions();

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckScalarResults(1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, kMaxIterations + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_diagonal, SolverTester,
                        *boost::unit_test::label("slow")) {
  dp3::ddecal::BdaIterativeDiagonalSolver solver;
  InitializeSolver(solver);
  solver.SetLLSSolverType(LLSSolverType::QR, 0.0, 0.0);

  BOOST_CHECK_EQUAL(solver.NSolutionPolarizations(), 2u);
  BOOST_REQUIRE_EQUAL(solver.ConstraintSolvers().size(), 1u);
  BOOST_CHECK_EQUAL(solver.ConstraintSolvers()[0], &solver);

  SetDiagonalSolutions();

  const dp3::ddecal::BdaSolverBuffer& solver_buffer = FillBDAData();
  dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks, kNDirections,
                              kNAntennas, Antennas1(), Antennas2());

  dp3::ddecal::SolverBase::SolveResult result =
      solver.Solve(data, GetSolverSolutions(), 0.0, nullptr);

  CheckDiagonalResults(1.0E-2);
  // The iterative solver solves the requested accuracy within the max
  // iterations so just check if the nr of iterations is <= max+1.
  BOOST_CHECK_LE(result.iterations, kMaxIterations + 1);
}

BOOST_AUTO_TEST_SUITE_END()
