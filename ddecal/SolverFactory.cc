// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverFactory.h"

#include "Settings.h"

#include "gain_solvers/BdaDiagonalSolver.h"
#include "gain_solvers/BdaHybridSolver.h"
#include "gain_solvers/BdaIterativeDiagonalSolver.h"
#include "gain_solvers/BdaIterativeScalarSolver.h"
#include "gain_solvers/BdaScalarSolver.h"

#include "gain_solvers/DiagonalSolver.h"
#include "gain_solvers/FullJonesSolver.h"
#include "gain_solvers/HybridSolver.h"
#include "gain_solvers/IterativeDiagonalSolver.h"
#include "gain_solvers/IterativeScalarSolver.h"
#include "gain_solvers/ScalarSolver.h"

#include "constraints/RotationConstraint.h"
#include "constraints/RotationAndDiagonalConstraint.h"
#ifdef HAVE_ARMADILLO
#include "constraints/ScreenConstraint.h"
#endif
#include "constraints/SmoothnessConstraint.h"
#include "constraints/TECConstraint.h"

#include <boost/make_unique.hpp>

namespace dp3 {
namespace ddecal {

namespace {

template <class SolverBaseType>
struct SolverTypes {};

template <>
struct SolverTypes<RegularSolverBase> {
  using Scalar = ScalarSolver;
  using IterativeScalar = IterativeScalarSolver;
  using Diagonal = DiagonalSolver;
  using IterativeDiagonal = IterativeDiagonalSolver;
  using Hybrid = HybridSolver;
};

template <>
struct SolverTypes<BdaSolverBase> {
  using Scalar = BdaScalarSolver;
  using IterativeScalar = BdaIterativeScalarSolver;
  using Diagonal = BdaDiagonalSolver;
  using IterativeDiagonal = BdaIterativeDiagonalSolver;
  using Hybrid = BdaHybridSolver;
};

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateScalarSolver(SolverAlgorithm algorithm) {
  if (ddecal::SolverAlgorithm::kDirectionIterative == algorithm)
    return boost::make_unique<
        typename SolverTypes<SolverBaseType>::IterativeScalar>();
  else
    return boost::make_unique<typename SolverTypes<SolverBaseType>::Scalar>();
}

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateDiagonalSolver(
    SolverAlgorithm algorithm) {
  if (ddecal::SolverAlgorithm::kDirectionIterative == algorithm)
    return boost::make_unique<
        typename SolverTypes<SolverBaseType>::IterativeDiagonal>();
  else
    return boost::make_unique<typename SolverTypes<SolverBaseType>::Diagonal>();
}

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateFullJonesSolver(
    SolverAlgorithm algorithm);

template <>
std::unique_ptr<RegularSolverBase> CreateFullJonesSolver(
    ddecal::SolverAlgorithm algorithm) {
  if (ddecal::SolverAlgorithm::kDirectionIterative == algorithm)
    throw std::runtime_error(
        "The direction-iterating algorithm is not available for the "
        "specified solving mode");
  else
    return boost::make_unique<FullJonesSolver>();
}

template <>
std::unique_ptr<BdaSolverBase> CreateFullJonesSolver([
    [maybe_unused]] ddecal::SolverAlgorithm algorithm) {
  throw std::runtime_error(
      "FullJones calibration not implemented in combination with BDA");
}

void AddConstraints(SolverBase& solver, const Settings& settings,
                    const common::ParameterSet& parset,
                    const std::string& prefix) {
  if (settings.core_constraint != 0.0 || !settings.antenna_constraint.empty()) {
    solver.AddConstraint(boost::make_unique<AntennaConstraint>());
  }
  if (settings.smoothness_constraint != 0.0) {
    solver.AddConstraint(boost::make_unique<SmoothnessConstraint>(
        settings.smoothness_constraint, settings.smoothness_ref_frequency));
  }

  switch (settings.mode) {
    case base::CalType::kScalar:
    case base::CalType::kDiagonal:
    case base::CalType::kFullJones:
      // no extra constraints
      break;
    case base::CalType::kScalarPhase:
    case base::CalType::kDiagonalPhase:
      solver.AddConstraint(boost::make_unique<ddecal::PhaseOnlyConstraint>());
      break;
    case base::CalType::kScalarAmplitude:
    case base::CalType::kDiagonalAmplitude:
      solver.AddConstraint(
          boost::make_unique<ddecal::AmplitudeOnlyConstraint>());
      break;
    case base::CalType::kTec:
    case base::CalType::kTecAndPhase: {
      const auto tec_mode = (settings.mode == base::CalType::kTec)
                                ? TECConstraint::TECOnlyMode
                                : TECConstraint::TECAndCommonScalarMode;
      std::unique_ptr<TECConstraint> constraint;

      if (settings.approximate_tec) {
        auto approxConstraint =
            boost::make_unique<ApproximateTECConstraint>(tec_mode);
        approxConstraint->SetMaxApproximatingIterations(
            settings.max_approx_iterations);
        approxConstraint->SetFittingChunkSize(settings.approx_chunk_size);
        constraint = std::move(approxConstraint);
      } else {
        constraint = boost::make_unique<TECConstraint>(tec_mode);
      }
      constraint->setDoPhaseReference(settings.phase_reference);
      solver.AddConstraint(std::move(constraint));
      break;
    }
#ifdef HAVE_ARMADILLO
    case base::CalType::kTecScreen:
      solver.AddConstraint(
          boost::make_unique<ScreenConstraint>(parset, prefix + "tecscreen."));
      break;
#endif
    case base::CalType::kRotationAndDiagonal: {
      auto constraint = boost::make_unique<RotationAndDiagonalConstraint>();
      constraint->SetDoRotationReference(settings.rotation_reference);
      solver.AddConstraint(std::move(constraint));
      break;
    }
    case base::CalType::kRotation:
      solver.AddConstraint(boost::make_unique<RotationConstraint>());
      break;
    default:
      throw std::runtime_error("Unexpected solving mode: " +
                               ToString(settings.mode));
  }
}

void InitializeSolver(SolverBase& solver, const Settings& settings) {
  solver.SetLLSSolverType(
      settings.lls_solver_type, settings.lls_min_tolerance,
      std::max(settings.lls_min_tolerance, settings.lls_max_tolerance));
  solver.SetMaxIterations(settings.max_iterations);
  solver.SetAccuracy(settings.tolerance);
  solver.SetConstraintAccuracy(settings.approx_tolerance);
  solver.SetStepSize(settings.step_size);
  solver.SetDetectStalling(settings.detect_stalling);
}

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateSolver(
    const Settings& settings, const common::ParameterSet& parset,
    const std::string& prefix, [[maybe_unused]] SolverAlgorithm algorithm) {
  std::unique_ptr<SolverBaseType> solver;
  switch (settings.mode) {
    case base::CalType::kScalar:
    case base::CalType::kScalarAmplitude:
      solver = CreateScalarSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kScalarPhase:
    case base::CalType::kTec:
    case base::CalType::kTecAndPhase:
      solver = CreateScalarSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(true);
      break;
    case base::CalType::kDiagonal:
    case base::CalType::kDiagonalAmplitude:
      solver = CreateDiagonalSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kDiagonalPhase:
      solver = CreateDiagonalSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(true);
      break;
    case base::CalType::kFullJones:
    case base::CalType::kRotationAndDiagonal:
    case base::CalType::kRotation:
      solver = CreateFullJonesSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kTecScreen:
#ifdef HAVE_ARMADILLO
      solver = CreateScalarSolver<SolverBaseType>(settings.solver_algorithm);
      solver->SetPhaseOnly(true);
#else
      throw std::runtime_error(
          "Can not use TEC screen: Armadillo is not available. Recompile DP3 "
          "with Armadillo.");
#endif
      break;
    default:
      throw std::runtime_error("Unexpected solving mode: " +
                               base::ToString(settings.mode));
  }

  AddConstraints(*solver, settings, parset, prefix);

  InitializeSolver(*solver, settings);

  return solver;
}

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateSolver(const Settings& settings,
                                             const common::ParameterSet& parset,
                                             const std::string& prefix) {
  std::unique_ptr<SolverBaseType> solver;

  if (settings.solver_algorithm == SolverAlgorithm::kHybrid) {
    auto a = CreateSolver<SolverBaseType>(settings, parset, prefix,
                                          SolverAlgorithm::kDirectionSolve);
    // The max_iterations is divided by 6 to use at most 1/6th of the iterations
    // in the first solver.
    a->SetMaxIterations(std::max<size_t>(1u, settings.max_iterations / 6u));
    auto b = CreateSolver<SolverBaseType>(settings, parset, prefix,
                                          SolverAlgorithm::kDirectionIterative);

    auto hybrid_solver =
        boost::make_unique<typename SolverTypes<SolverBaseType>::Hybrid>();
    hybrid_solver->SetMaxIterations(settings.max_iterations);
    hybrid_solver->AddSolver(std::move(a));
    hybrid_solver->AddSolver(std::move(b));
    solver = std::move(hybrid_solver);
  } else {
    solver = CreateSolver<SolverBaseType>(settings, parset, prefix,
                                          settings.solver_algorithm);
  }
  return solver;
}

}  // namespace

std::unique_ptr<RegularSolverBase> CreateRegularSolver(
    const Settings& settings, const common::ParameterSet& parset,
    const std::string& prefix) {
  return CreateSolver<RegularSolverBase>(settings, parset, prefix);
}

std::unique_ptr<BdaSolverBase> CreateBdaSolver(
    const Settings& settings, const common::ParameterSet& parset,
    const std::string& prefix) {
  return CreateSolver<BdaSolverBase>(settings, parset, prefix);
}

}  // namespace ddecal
}  // namespace dp3