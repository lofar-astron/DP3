// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SolverFactory.h"

#include "../../Settings.h"
#include "../../constraints/RotationAndDiagonalConstraint.h"
#include "../../constraints/RotationConstraint.h"

#ifdef HAVE_ARMADILLO
#include "../../constraints/ScreenConstraint.h"
#endif

#include "../../constraints/SmoothnessConstraint.h"
#include "../../constraints/TECConstraint.h"
#include "../../gain_solvers/BdaDiagonalSolver.h"
#include "../../gain_solvers/BdaIterativeScalarSolver.h"
#include "../../gain_solvers/BdaIterativeDiagonalSolver.h"
#include "../../gain_solvers/BdaScalarSolver.h"
#include "../../gain_solvers/DiagonalSolver.h"
#include "../../gain_solvers/FullJonesSolver.h"
#include "../../gain_solvers/IterativeScalarSolver.h"
#include "../../gain_solvers/IterativeDiagonalSolver.h"
#include "../../gain_solvers/ScalarSolver.h"

#include "../../../common/ParameterSet.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::CreateBdaSolver;
using dp3::ddecal::CreateRegularSolver;
using dp3::ddecal::Settings;

namespace {

dp3::common::ParameterSet ParsetForMode(const std::string& mode) {
  dp3::common::ParameterSet parset;
  parset.add("msin", "foo");  // msin is a mandatory field for Settings.
  parset.add("mode", mode);
  return parset;
}

template <class ExpectedType>
void CheckRegularSolverType(const dp3::common::ParameterSet& parset) {
  const Settings settings(parset, "");
  std::unique_ptr<dp3::ddecal::RegularSolverBase> solver =
      CreateRegularSolver(settings, parset, "");
  BOOST_CHECK(dynamic_cast<ExpectedType*>(solver.get()));
}

template <class ExpectedType>
void CheckBdaSolverType(const dp3::common::ParameterSet& parset) {
  const Settings settings(parset, "");
  std::unique_ptr<dp3::ddecal::BdaSolverBase> solver =
      CreateBdaSolver(settings, parset, "");
  BOOST_CHECK(dynamic_cast<ExpectedType*>(solver.get()));
}

template <class ExpectedType>
void CheckConstraintType(dp3::ddecal::SolverBase& solver) {
  const std::vector<std::unique_ptr<dp3::ddecal::Constraint>>& constraints =
      solver.GetConstraints();
  BOOST_REQUIRE_EQUAL(constraints.size(), 1u);
  BOOST_CHECK(dynamic_cast<ExpectedType*>(constraints.front().get()));
}

template <class ExpectedType>
void CheckConstraintType(const dp3::common::ParameterSet& parset) {
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  CheckConstraintType<ExpectedType>(*solver);

  auto bda_solver = CreateBdaSolver(settings, parset, "");
  CheckConstraintType<ExpectedType>(*bda_solver);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(solverfactory)

#ifndef HAVE_ARMADILLO
BOOST_DATA_TEST_CASE(scalar_type,
                     boost::unit_test::data::make({"scalar", "scalaramplitude",
                                                   "scalarphase", "tec",
                                                   "tecandphase"}),
                     mode) {
#else
BOOST_DATA_TEST_CASE(scalar_type,
                     boost::unit_test::data::make({"scalar", "scalaramplitude",
                                                   "scalarphase", "tec",
                                                   "tecandphase", "tecscreen"}),
                     mode) {
#endif
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  CheckRegularSolverType<dp3::ddecal::ScalarSolver>(parset);
  CheckBdaSolverType<dp3::ddecal::BdaScalarSolver>(parset);

  parset.add("solveralgorithm", "directioniterative");
  CheckRegularSolverType<dp3::ddecal::IterativeScalarSolver>(parset);
  CheckBdaSolverType<dp3::ddecal::BdaIterativeScalarSolver>(parset);
}

BOOST_DATA_TEST_CASE(diagonal_type,
                     boost::unit_test::data::make(
                         {"diagonal", "diagonalamplitude", "diagonalphase"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  CheckRegularSolverType<dp3::ddecal::DiagonalSolver>(parset);
  CheckBdaSolverType<dp3::ddecal::BdaDiagonalSolver>(parset);

  parset.add("solveralgorithm", "directioniterative");
  CheckRegularSolverType<dp3::ddecal::IterativeDiagonalSolver>(parset);
  CheckBdaSolverType<dp3::ddecal::BdaIterativeDiagonalSolver>(parset);
}

BOOST_DATA_TEST_CASE(fulljones_type,
                     boost::unit_test::data::make(
                         {"fulljones", "rotation+diagonal", "rotation"}),
                     mode) {
  // Currently, only the regular non-iterative full-jones solver is implemented.
  // Requesting other full jones solvers (bda / iterative) raises an exception.
  // Please also update other tests below when more full jones solves become
  // available, since they now skip those tests.

  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings_plain(parset, "");
  CheckRegularSolverType<dp3::ddecal::FullJonesSolver>(parset);

  BOOST_CHECK_THROW(dp3::ddecal::CreateBdaSolver(settings_plain, parset, ""),
                    std::runtime_error);

  parset.add("solveralgorithm", "directioniterative");
  const Settings settings_iterative(parset, "");
  BOOST_CHECK_THROW(CreateRegularSolver(settings_iterative, parset, ""),
                    std::runtime_error);
  BOOST_CHECK_THROW(CreateBdaSolver(settings_iterative, parset, ""),
                    std::runtime_error);
}

BOOST_DATA_TEST_CASE(
    no_phase_only,
    boost::unit_test::data::make({"scalar", "scalaramplitude", "diagonal",
                                  "diagonalamplitude", "fulljones",
                                  "rotation+diagonal", "rotation"}),
    mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  BOOST_CHECK(!solver->GetPhaseOnly());
}

BOOST_DATA_TEST_CASE(no_phase_only_bda,
                     boost::unit_test::data::make({"scalar", "scalaramplitude",
                                                   "diagonal",
                                                   "diagonalamplitude"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");

  auto bda_solver = CreateBdaSolver(settings, parset, "");
  BOOST_CHECK(!bda_solver->GetPhaseOnly());
}

#ifndef HAVE_ARMADILLO
BOOST_DATA_TEST_CASE(phase_only,
                     boost::unit_test::data::make({"scalarphase",
                                                   "diagonalphase", "tec",
                                                   "tecandphase"}),
                     mode) {
#else
BOOST_DATA_TEST_CASE(phase_only,
                     boost::unit_test::data::make({"scalarphase",
                                                   "diagonalphase", "tec",
                                                   "tecandphase", "tecscreen"}),
                     mode) {
#endif
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  BOOST_CHECK(solver->GetPhaseOnly());

  auto bda_solver = CreateBdaSolver(settings, parset, "");
  BOOST_CHECK(bda_solver->GetPhaseOnly());
}

BOOST_DATA_TEST_CASE(no_constraints,
                     boost::unit_test::data::make({"scalar", "diagonal",
                                                   "fulljones"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");
  auto solver = CreateRegularSolver(settings, parset, "");
  BOOST_CHECK(solver->GetConstraints().empty());

  if (mode != std::string("fulljones")) {
    auto bda_solver = CreateBdaSolver(settings, parset, "");
    BOOST_CHECK(bda_solver->GetConstraints().empty());
  }
}

BOOST_DATA_TEST_CASE(phase_constraint,
                     boost::unit_test::data::make({"scalarphase",
                                                   "diagonalphase"}),
                     mode) {
  CheckConstraintType<dp3::ddecal::PhaseOnlyConstraint>(ParsetForMode(mode));
}

BOOST_DATA_TEST_CASE(amplitude_constraint,
                     boost::unit_test::data::make({"scalaramplitude",
                                                   "diagonalamplitude"}),
                     mode) {
  CheckConstraintType<dp3::ddecal::AmplitudeOnlyConstraint>(
      ParsetForMode(mode));
}

BOOST_DATA_TEST_CASE(tec_constraint,
                     boost::unit_test::data::make({"tec", "tecandphase"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings_default(parset, "");

  // For the regular solver, use default settings.
  auto solver = CreateRegularSolver(settings_default, parset, "");
  CheckConstraintType<dp3::ddecal::TECConstraint>(*solver);
  BOOST_CHECK(!dynamic_cast<dp3::ddecal::ApproximateTECConstraint*>(
      solver->GetConstraints()[0].get()));

  // For the BDA solver, use an approximate tec constraint.
  parset.add("approximatetec", "true");
  const Settings settings_approx(parset, "");

  auto bda_solver = CreateBdaSolver(settings_approx, parset, "");
  CheckConstraintType<dp3::ddecal::ApproximateTECConstraint>(*bda_solver);
}

#ifdef HAVE_ARMADILLO
BOOST_AUTO_TEST_CASE(screen_constraint) {
  CheckConstraintType<dp3::ddecal::ScreenConstraint>(
      ParsetForMode("tecscreen"));
}
#endif

BOOST_AUTO_TEST_CASE(rotation_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("rotation");
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  CheckConstraintType<dp3::ddecal::RotationConstraint>(*solver);
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("rotation+diagonal");
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  CheckConstraintType<dp3::ddecal::RotationAndDiagonalConstraint>(*solver);
}

BOOST_AUTO_TEST_CASE(core_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("diagonal");
  parset.add("coreconstraint", "0.42");
  CheckConstraintType<dp3::ddecal::AntennaConstraint>(parset);
}

BOOST_AUTO_TEST_CASE(antenna_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("scalar");
  parset.add("antennaconstraint", "[[foo,bar]]");
  CheckConstraintType<dp3::ddecal::AntennaConstraint>(parset);
}

BOOST_AUTO_TEST_CASE(smoothness_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("scalar");
  parset.add("smoothnessconstraint", "0.42");
  CheckConstraintType<dp3::ddecal::SmoothnessConstraint>(parset);
}

BOOST_AUTO_TEST_CASE(multiple_constraints) {
  dp3::common::ParameterSet parset = ParsetForMode("diagonalphase");
  parset.add("antennaconstraint", "[[foo,bar]]");
  parset.add("smoothnessconstraint", "0.42");
  const Settings settings(parset, "");

  auto solver = CreateRegularSolver(settings, parset, "");
  const std::vector<std::unique_ptr<dp3::ddecal::Constraint>>& constraints =
      solver->GetConstraints();
  BOOST_REQUIRE_EQUAL(constraints.size(), 3u);
  BOOST_CHECK(
      dynamic_cast<dp3::ddecal::AntennaConstraint*>(constraints[0].get()));
  BOOST_CHECK(
      dynamic_cast<dp3::ddecal::SmoothnessConstraint*>(constraints[1].get()));
  BOOST_CHECK(
      dynamic_cast<dp3::ddecal::PhaseOnlyConstraint*>(constraints[2].get()));
}

BOOST_AUTO_TEST_CASE(solver_settings_default) {
  dp3::common::ParameterSet parset = ParsetForMode("tec");
  parset.add("approximatetec", "true");
  const Settings settings(parset, "");
  auto solver = CreateRegularSolver(settings, parset, "");

  BOOST_CHECK(solver->GetLLSSolverType() == dp3::ddecal::LLSSolverType::QR);
  std::pair<double, double> tolerance = solver->GetLLSSolverTolerance();
  BOOST_CHECK_CLOSE(tolerance.first, 1.0e-7, 1.0e-8);
  BOOST_CHECK_CLOSE(tolerance.second, 1.0e-7, 1.0e-8);

  BOOST_CHECK_EQUAL(solver->GetMaxIterations(), 50u);
  BOOST_CHECK_CLOSE(solver->GetAccuracy(), 1.0e-4, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetConstraintAccuracy(), 1.0e-3, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetStepSize(), 0.2, 1.0e-8);
  BOOST_CHECK_EQUAL(solver->GetDetectStalling(), true);
}

BOOST_AUTO_TEST_CASE(solver_settings_custom) {
  dp3::common::ParameterSet parset = ParsetForMode("tec");
  parset.add("approximatetec", "true");

  parset.add("llssolver", "lsmr");
  parset.add("llsstarttolerance", "41e-6");
  parset.add("llstolerance", "42e-6");
  parset.add("maxiter", "42");
  parset.add("tolerance", "0.042");
  parset.add("approxtolerance", "0.1337");
  parset.add("stepsize", "0.42");
  parset.add("detectstalling", "false");
  const Settings settings(parset, "");
  auto solver = CreateRegularSolver(settings, parset, "");

  BOOST_CHECK(solver->GetLLSSolverType() == dp3::ddecal::LLSSolverType::LSMR);
  std::pair<double, double> tolerance = solver->GetLLSSolverTolerance();
  BOOST_CHECK_CLOSE(tolerance.first, 41e-6, 1.0e-8);
  BOOST_CHECK_CLOSE(tolerance.second, 42e-6, 1.0e-8);

  BOOST_CHECK_EQUAL(solver->GetMaxIterations(), 42u);
  BOOST_CHECK_CLOSE(solver->GetAccuracy(), 0.042, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetConstraintAccuracy(), 0.1337, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetStepSize(), 0.42, 1.0e-8);
  BOOST_CHECK_EQUAL(solver->GetDetectStalling(), false);
}

BOOST_AUTO_TEST_SUITE_END()
