// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SolverFactory.h"

#include "../../Settings.h"

#include "../../constraints/AntennaConstraint.h"
#include "../../constraints/AmplitudeOnlyConstraint.h"
#include "../../constraints/PhaseOnlyConstraint.h"
#include "../../constraints/RotationAndDiagonalConstraint.h"
#include "../../constraints/RotationConstraint.h"

#ifdef HAVE_ARMADILLO
#include "../../constraints/ScreenConstraint.h"
#endif

#include "../../constraints/SmoothnessConstraint.h"
#include "../../constraints/TECConstraint.h"
#include "../../gain_solvers/DiagonalSolver.h"
#include "../../gain_solvers/FullJonesSolver.h"
#include "../../gain_solvers/IterativeDiagonalSolver.h"
#include "../../gain_solvers/IterativeFullJonesSolver.h"
#include "../../gain_solvers/IterativeScalarSolver.h"
#include "../../gain_solvers/ScalarSolver.h"

#include "../../../common/ParameterSet.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::AntennaConstraint;
using dp3::ddecal::Constraint;
using dp3::ddecal::CreateSolver;
using dp3::ddecal::InitializeSolverConstraints;
using dp3::ddecal::PhaseOnlyConstraint;
using dp3::ddecal::Settings;
using dp3::ddecal::SmoothnessConstraint;
using dp3::ddecal::SolverBase;

namespace {

dp3::common::ParameterSet ParsetForMode(const std::string& mode) {
  dp3::common::ParameterSet parset;
  parset.add("msin", "foo");  // msin is a mandatory field for Settings.
  parset.add("mode", mode);
  return parset;
}

template <class ExpectedType>
void CheckSolverType(const dp3::common::ParameterSet& parset) {
  const Settings settings(parset, "");
  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  BOOST_CHECK(dynamic_cast<ExpectedType*>(solver.get()));
}

template <class ExpectedType>
const ExpectedType& CheckConstraintType(SolverBase& solver) {
  const std::vector<std::unique_ptr<Constraint>>& constraints =
      solver.GetConstraints();
  BOOST_REQUIRE_EQUAL(constraints.size(), 1u);
  ExpectedType* constraint =
      dynamic_cast<ExpectedType*>(constraints.front().get());
  BOOST_CHECK(constraint);
  return *constraint;
}

template <class ExpectedType>
void CheckConstraintType(const dp3::common::ParameterSet& parset) {
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  CheckConstraintType<ExpectedType>(*solver);
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
  CheckSolverType<dp3::ddecal::ScalarSolver>(parset);

  parset.add("solveralgorithm", "directioniterative");
  CheckSolverType<dp3::ddecal::IterativeScalarSolver>(parset);
}

BOOST_DATA_TEST_CASE(diagonal_type,
                     boost::unit_test::data::make(
                         {"diagonal", "diagonalamplitude", "diagonalphase"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  CheckSolverType<dp3::ddecal::DiagonalSolver>(parset);

  parset.add("solveralgorithm", "directioniterative");
  CheckSolverType<dp3::ddecal::IterativeDiagonalSolver>(parset);
}

BOOST_DATA_TEST_CASE(full_jones_type,
                     boost::unit_test::data::make(
                         {"fulljones", "rotation+diagonal", "rotation"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  CheckSolverType<dp3::ddecal::FullJonesSolver>(parset);

  parset.add("solveralgorithm", "directioniterative");
  CheckSolverType<dp3::ddecal::IterativeFullJonesSolver>(parset);
}

BOOST_DATA_TEST_CASE(
    no_phase_only,
    boost::unit_test::data::make({"scalar", "scalaramplitude", "diagonal",
                                  "diagonalamplitude", "fulljones",
                                  "rotation+diagonal", "rotation"}),
    mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  BOOST_CHECK(!solver->GetPhaseOnly());
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

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  BOOST_CHECK(solver->GetPhaseOnly());
}

BOOST_DATA_TEST_CASE(no_constraints,
                     boost::unit_test::data::make({"scalar", "diagonal",
                                                   "fulljones"}),
                     mode) {
  dp3::common::ParameterSet parset = ParsetForMode(mode);
  const Settings settings(parset, "");
  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  BOOST_CHECK(solver->GetConstraints().empty());
}

BOOST_DATA_TEST_CASE(phase_constraint,
                     boost::unit_test::data::make({"scalarphase",
                                                   "diagonalphase"}),
                     mode) {
  CheckConstraintType<PhaseOnlyConstraint>(ParsetForMode(mode));
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
  // For this test, only the size of these vectors matters.
  const std::vector<std::array<double, 3>> kAntennaPositions(5);
  const std::vector<std::string> kAntennaNames(kAntennaPositions.size());
  const size_t kNDirections = 4;
  const std::vector<dp3::base::Direction> kSourceDirections(kNDirections);
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);
  const std::vector<double> kFrequencies(10);

  dp3::common::ParameterSet parset = ParsetForMode(mode);

  {
    // Test a solver with default settings.
    const Settings settings(parset, "");
    std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
    InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                                kAntennaNames, kSolutionsPerDirecton,
                                kSourceDirections, kFrequencies);

    const Constraint& constraint =
        CheckConstraintType<dp3::ddecal::TECConstraint>(*solver);
    BOOST_CHECK(!dynamic_cast<dp3::ddecal::ApproximateTECConstraint*>(
        solver->GetConstraints()[0].get()));
    BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
    BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
    BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  }

  {
    // Test a solver with custom settings.
    parset.add("approximatetec", "true");
    const Settings settings(parset, "");
    std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
    InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                                kAntennaNames, kSolutionsPerDirecton,
                                kSourceDirections, kFrequencies);

    const Constraint& constraint =
        CheckConstraintType<dp3::ddecal::ApproximateTECConstraint>(*solver);
    BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
    BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
    BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  }
}

#ifdef HAVE_ARMADILLO
BOOST_AUTO_TEST_CASE(screen_constraint) {
  using dp3::ddecal::PiercePoint;
  using dp3::ddecal::ScreenConstraint;

  const std::vector<std::array<double, 3>> kAntennaPositions{
      {1.0, 1.0, 1.0}, {100.0, 1.0, 1.0}, {1.0, 1.0, 2.0}};
  const std::vector<std::string> kAntennaNames(3);
  const size_t kNDirections = 2;
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);
  const std::vector<dp3::base::Direction> kSourceDirections{{0.1, 0.1},
                                                            {0.3, 0.3}};
  const std::vector<double> kFrequencies{42.0, 43.0, 44.0, 45.0};
  const double kHeight = 42.0e3;
  const std::vector<std::size_t> kExpectedCoreAntennas{0, 2};

  dp3::common::ParameterSet parset = ParsetForMode("tecscreen");
  parset.add("tecscreen.coreconstraint", "42");
  parset.add("tecscreen.height", std::to_string(kHeight));
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                              kAntennaNames, kSolutionsPerDirecton,
                              kSourceDirections, kFrequencies);

  const ScreenConstraint& constraint =
      CheckConstraintType<ScreenConstraint>(*solver);
  BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
  BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
  BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      constraint.GetCoreAntennas().begin(), constraint.GetCoreAntennas().end(),
      kExpectedCoreAntennas.begin(), kExpectedCoreAntennas.end());
  const std::vector<std::vector<PiercePoint>>& pierce_points =
      constraint.GetPiercePoints();
  BOOST_REQUIRE_EQUAL(pierce_points.size(), kAntennaPositions.size());
  for (size_t ant = 0; ant < kAntennaPositions.size(); ++ant) {
    BOOST_REQUIRE_EQUAL(pierce_points[ant].size(), kSourceDirections.size());
    for (size_t dir = 0; dir < kSourceDirections.size(); ++dir) {
      const PiercePoint& pp = pierce_points[ant][dir];
      BOOST_CHECK_EQUAL(pp.getPos().getValue(),
                        casacore::MVPosition(kAntennaPositions[ant][0],
                                             kAntennaPositions[ant][1],
                                             kAntennaPositions[ant][2]));
      BOOST_CHECK_EQUAL(pp.getDir().getValue(),
                        casacore::MVDirection(kSourceDirections[dir].ra,
                                              kSourceDirections[dir].dec));
      BOOST_CHECK_EQUAL(pp.getHeight(), kHeight);
    }
  }
}
#endif

BOOST_AUTO_TEST_CASE(rotation_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("rotation");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  CheckConstraintType<dp3::ddecal::RotationConstraint>(*solver);
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal_constraint) {
  dp3::common::ParameterSet parset = ParsetForMode("rotation+diagonal");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  CheckConstraintType<dp3::ddecal::RotationAndDiagonalConstraint>(*solver);
}

BOOST_AUTO_TEST_CASE(core_constraint) {
  // Only antennas 0 and 2 are within the coreconstraint.
  const std::vector<std::array<double, 3>> kAntennaPositions{
      {1.0, 1.0, 1.0}, {100.0, 1.0, 1.0}, {1.0, 1.0, 2.0}};
  const std::vector<std::set<size_t>> kExpectedAntennaSets{{0, 2}};
  const std::vector<std::string> kAntennaNames(3);
  const size_t kNDirections = 5;
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);
  const std::vector<dp3::base::Direction> kSourceDirections(kNDirections);
  const std::vector<double> kFrequencies{42.0, 43.0, 44.0, 45.0};

  dp3::common::ParameterSet parset = ParsetForMode("diagonal");
  parset.add("coreconstraint", "42");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                              kAntennaNames, kSolutionsPerDirecton,
                              kSourceDirections, kFrequencies);

  const AntennaConstraint& constraint =
      CheckConstraintType<AntennaConstraint>(*solver);
  BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
  BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
  BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  BOOST_CHECK(constraint.GetAntennaSets() == kExpectedAntennaSets);
}

BOOST_AUTO_TEST_CASE(antenna_constraint) {
  const std::vector<std::array<double, 3>> kAntennaPositions(6);
  const std::vector<std::string> kAntennaNames{"Jan", "foo",  "X",
                                               "bar", "Piet", "7TiMeS6"};
  const std::vector<std::set<size_t>> kExpectedAntennaSets{
      {1, 3, 5}, {1, 4}, {0, 3}};
  const size_t kNDirections = 5;
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);
  const std::vector<dp3::base::Direction> kSourceDirections(kNDirections);
  const std::vector<double> kFrequencies{42.0};

  dp3::common::ParameterSet parset = ParsetForMode("scalar");
  parset.add("antennaconstraint", "[[foo,bar,7TiMeS6],[Piet,foo],[bar,Jan]]");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                              kAntennaNames, kSolutionsPerDirecton,
                              kSourceDirections, kFrequencies);

  const AntennaConstraint& constraint =
      CheckConstraintType<AntennaConstraint>(*solver);
  BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
  BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
  BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  BOOST_CHECK(constraint.GetAntennaSets() == kExpectedAntennaSets);
}

BOOST_AUTO_TEST_CASE(antenna_constraint_invalid) {
  const std::vector<std::array<double, 3>> kAntennaPositions(2);
  const std::vector<std::string> kAntennaNames{"foo", "bar"};
  const size_t kNDirections = 5;
  const std::vector<dp3::base::Direction> kSourceDirections(kNDirections);
  const std::vector<double> kFrequencies{42.0};
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);

  // Settings::ReadAntennaConstraint already catches syntax errors in the parset
  // value. The SolverFactory catches semantic errors, such as invalid names.
  dp3::common::ParameterSet parset = ParsetForMode("scalar");
  parset.add("antennaconstraint", "[[foo, bar, missing]]");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  CheckConstraintType<AntennaConstraint>(*solver);

  BOOST_CHECK_THROW(InitializeSolverConstraints(
                        *solver, settings, kAntennaPositions, kAntennaNames,
                        kSolutionsPerDirecton, kSourceDirections, kFrequencies),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(smoothness_constraint_without_ref_distance) {
  const std::vector<std::array<double, 3>> kAntennaPositions{
      {1, 1, 1}, {2, 1, 1}, {1, 3, 1}, {1, 1, 4}};
  const std::vector<std::string> kAntennaNames(kAntennaPositions.size());
  const size_t kNDirections = 1;
  const std::vector<dp3::base::Direction> kSourceDirections(1);
  const std::vector<double> kFrequencies{42.0};
  const std::vector<double> kExpectedDistanceFactors(kAntennaPositions.size(),
                                                     1.0);
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);

  dp3::common::ParameterSet parset = ParsetForMode("scalar");
  parset.add("smoothnessconstraint", "1");
  parset.add("smoothnessreffrequency", "200");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                              kAntennaNames, kSolutionsPerDirecton,
                              kSourceDirections, kFrequencies);

  const SmoothnessConstraint& constraint =
      CheckConstraintType<SmoothnessConstraint>(*solver);
  BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
  BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
  BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  BOOST_CHECK(constraint.GetDistanceFactors() == kExpectedDistanceFactors);
}

BOOST_AUTO_TEST_CASE(smoothness_constraint_with_ref_distance) {
  const std::vector<std::array<double, 3>> kAntennaPositions{
      {1, 1, 1}, {2, 1, 1}, {1, 7, 1}, {1, 1, 8}};
  const std::vector<std::string> kAntennaNames(kAntennaPositions.size());
  const std::vector<dp3::base::Direction> kSourceDirections(1);
  const size_t kNDirections = 1;
  const std::vector<double> kFrequencies{142.0};
  const std::vector<double> kExpectedDistanceFactors{42, 42, 7, 6};
  const std::vector<size_t> kSolutionsPerDirecton(kNDirections, 1);

  dp3::common::ParameterSet parset = ParsetForMode("diagonal");
  parset.add("smoothnessconstraint", "1");
  parset.add("smoothnessreffrequency", "200");
  parset.add("smoothnessrefdistance", "42");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  InitializeSolverConstraints(*solver, settings, kAntennaPositions,
                              kAntennaNames, kSolutionsPerDirecton,
                              kSourceDirections, kFrequencies);

  const SmoothnessConstraint& constraint =
      CheckConstraintType<SmoothnessConstraint>(*solver);
  BOOST_CHECK_EQUAL(constraint.NAntennas(), kAntennaNames.size());
  BOOST_CHECK_EQUAL(constraint.NDirections(), kSourceDirections.size());
  BOOST_CHECK_EQUAL(constraint.NChannelBlocks(), kFrequencies.size());
  BOOST_CHECK(constraint.GetDistanceFactors() == kExpectedDistanceFactors);
}

BOOST_AUTO_TEST_CASE(multiple_constraints) {
  dp3::common::ParameterSet parset = ParsetForMode("diagonalphase");
  parset.add("antennaconstraint", "[[foo,bar]]");
  parset.add("smoothnessconstraint", "0.42");
  const Settings settings(parset, "");

  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");
  const std::vector<std::unique_ptr<dp3::ddecal::Constraint>>& constraints =
      solver->GetConstraints();
  BOOST_REQUIRE_EQUAL(constraints.size(), 3u);
  BOOST_CHECK(dynamic_cast<AntennaConstraint*>(constraints[0].get()));
  BOOST_CHECK(dynamic_cast<SmoothnessConstraint*>(constraints[1].get()));
  BOOST_CHECK(dynamic_cast<PhaseOnlyConstraint*>(constraints[2].get()));
}

BOOST_AUTO_TEST_CASE(solver_settings_default) {
  dp3::common::ParameterSet parset = ParsetForMode("tec");
  parset.add("approximatetec", "true");
  const Settings settings(parset, "");
  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");

  BOOST_CHECK(solver->GetLLSSolverType() == dp3::ddecal::LLSSolverType::QR);

  BOOST_CHECK_EQUAL(solver->GetMaxIterations(), 50u);
  BOOST_CHECK_CLOSE(solver->GetAccuracy(), 1.0e-4, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetConstraintAccuracy(), 1.0e-3, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetStepSize(), 0.2, 1.0e-8);
  BOOST_CHECK_EQUAL(solver->GetDetectStalling(), true);
}

BOOST_AUTO_TEST_CASE(solver_settings_custom) {
  dp3::common::ParameterSet parset = ParsetForMode("tec");
  parset.add("approximatetec", "true");

  parset.add("llssolver", "normalequations");
  parset.add("llsstarttolerance", "41e-6");
  parset.add("llstolerance", "42e-6");
  parset.add("maxiter", "42");
  parset.add("tolerance", "0.042");
  parset.add("approxtolerance", "0.1337");
  parset.add("stepsize", "0.42");
  parset.add("detectstalling", "false");
  const Settings settings(parset, "");
  std::unique_ptr<SolverBase> solver = CreateSolver(settings, parset, "");

  BOOST_CHECK(solver->GetLLSSolverType() ==
              dp3::ddecal::LLSSolverType::NORMAL_EQUATIONS);

  BOOST_CHECK_EQUAL(solver->GetMaxIterations(), 42u);
  BOOST_CHECK_CLOSE(solver->GetAccuracy(), 0.042, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetConstraintAccuracy(), 0.1337, 1.0e-8);
  BOOST_CHECK_CLOSE(solver->GetStepSize(), 0.42, 1.0e-8);
  BOOST_CHECK_EQUAL(solver->GetDetectStalling(), false);
}

BOOST_AUTO_TEST_SUITE_END()
