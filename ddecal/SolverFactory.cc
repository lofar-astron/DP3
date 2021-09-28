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
  if (SolverAlgorithm::kDirectionIterative == algorithm)
    return boost::make_unique<
        typename SolverTypes<SolverBaseType>::IterativeScalar>();
  else
    return boost::make_unique<typename SolverTypes<SolverBaseType>::Scalar>();
}

template <class SolverBaseType>
std::unique_ptr<SolverBaseType> CreateDiagonalSolver(
    SolverAlgorithm algorithm) {
  if (SolverAlgorithm::kDirectionIterative == algorithm)
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
    SolverAlgorithm algorithm) {
  if (SolverAlgorithm::kDirectionIterative == algorithm)
    throw std::runtime_error(
        "The direction-iterating algorithm is not available for the "
        "specified solving mode");
  else
    return boost::make_unique<FullJonesSolver>();
}

template <>
std::unique_ptr<BdaSolverBase> CreateFullJonesSolver(
    [[maybe_unused]] SolverAlgorithm algorithm) {
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
      solver.AddConstraint(boost::make_unique<PhaseOnlyConstraint>());
      break;
    case base::CalType::kScalarAmplitude:
    case base::CalType::kDiagonalAmplitude:
      solver.AddConstraint(boost::make_unique<AmplitudeOnlyConstraint>());
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

/**
 * Find all antennas / stations within certain distance.
 * @param core_constraint Maximum distance for core antennas.
 * @param positions Positions of the antennas. The function uses the
 * the first antenna as reference antenna.
 * @return The indices of the antennas that are within the maximum distance
 * of the reference station.
 */
std::set<size_t> DetermineCoreAntennas(
    double core_constraint,
    const std::vector<std::array<double, 3>>& positions) {
  std::set<size_t> core_indices;
  if (!positions.empty()) {
    const double refx = positions[0][0];
    const double refy = positions[0][1];
    const double refz = positions[0][2];
    const double core_dist_squared = core_constraint * core_constraint;
    for (size_t ant = 0; ant != positions.size(); ++ant) {
      const double dx = refx - positions[ant][0];
      const double dy = refy - positions[ant][1];
      const double dz = refz - positions[ant][2];
      const double dist_squared = dx * dx + dy * dy + dz * dz;
      if (dist_squared <= core_dist_squared) {
        core_indices.insert(core_indices.end(), ant);
      }
    }
  }
  return core_indices;
}

void InitializeAntennaCoreConstraint(
    AntennaConstraint& constraint, double core_constraint,
    const std::vector<std::array<double, 3>>& antenna_positions) {
  std::vector<std::set<size_t>> constraint_list;
  constraint_list.push_back(
      DetermineCoreAntennas(core_constraint, antenna_positions));
  constraint.SetAntennaSets(std::move(constraint_list));
}

void InitializeAntennaConstraint(
    AntennaConstraint& constraint,
    const std::vector<std::set<std::string>>& constraint_name_groups,
    const std::vector<std::string>& antenna_names) {
  // Set the antenna constraint to a list of stations indices that
  // are to be kept the same during the solve.
  std::vector<std::set<size_t>> constraint_list;
  for (const std::set<std::string>& constraint_name_set :
       constraint_name_groups) {
    if (constraint_name_set.size() <= 1) {
      throw std::runtime_error(
          "Error in antenna constraint: at least two antennas expected");
    }

    constraint_list.emplace_back();
    for (const std::string& constraint_name : constraint_name_set) {
      const auto iter = std::find(antenna_names.begin(), antenna_names.end(),
                                  constraint_name);
      if (iter == antenna_names.end()) {
        throw std::runtime_error("Error in antenna constraint: Antenna '" +
                                 constraint_name + "' not found.");
      }
      constraint_list.back().insert(constraint_list.back().end(),
                                    iter - antenna_names.begin());
    }
  }
  constraint.SetAntennaSets(std::move(constraint_list));
}

#ifdef HAVE_ARMADILLO
void InitializeScreenConstraint(
    ScreenConstraint& constraint, double core_constraint,
    const std::vector<std::array<double, 3>>& antenna_positions,
    const std::vector<std::pair<double, double>>& source_positions) {
  constraint.InitPiercePoints(antenna_positions, source_positions);
  constraint.SetCoreAntennas(
      DetermineCoreAntennas(core_constraint, antenna_positions));
}
#endif

void InitializeSmoothnessConstraint(
    SmoothnessConstraint& constraint, double ref_distance,
    const std::vector<std::array<double, 3>>& antenna_positions) {
  std::vector<double> distance_factors;
  // If no smoothness reference distance is specified, the smoothing is
  // made independent of the distance
  if (ref_distance == 0.0) {
    distance_factors.assign(antenna_positions.size(), 1.0);
  } else {
    // Make a list of factors such that more distant antennas apply a
    // smaller smoothing kernel.
    distance_factors.reserve(antenna_positions.size());
    for (size_t i = 1; i != antenna_positions.size(); ++i) {
      const double dx = antenna_positions[0][0] - antenna_positions[i][0];
      const double dy = antenna_positions[0][1] - antenna_positions[i][1];
      const double dz = antenna_positions[0][2] - antenna_positions[i][2];
      const double factor =
          ref_distance / std::sqrt(dx * dx + dy * dy + dz * dz);
      distance_factors.push_back(factor);
      // For antenna 0, the distance of antenna 1 is used:
      if (i == 1) distance_factors.push_back(factor);
    }
  }
  constraint.SetDistanceFactors(std::move(distance_factors));
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

void InitializeSolverConstraints(
    SolverBase& solver, const Settings& settings,
    const std::vector<std::array<double, 3>>& antenna_positions,
    const std::vector<std::string>& antenna_names,
    const std::vector<std::pair<double, double>>& source_positions,
    const std::vector<double>& frequencies) {
  for (const std::unique_ptr<Constraint>& constraint :
       solver.GetConstraints()) {
    // Initialize the constraint with some common metadata.
    constraint->Initialize(antenna_positions.size(), source_positions.size(),
                           frequencies);

    // Different constraints need different information. Determine if the
    // constraint is of a type that needs more information, and if so
    // initialize the constraint.
    AntennaConstraint* antenna_constraint =
        dynamic_cast<AntennaConstraint*>(constraint.get());
    if (antenna_constraint) {
      if (settings.antenna_constraint.empty()) {
        InitializeAntennaCoreConstraint(
            *antenna_constraint, settings.core_constraint, antenna_positions);
      } else {
        InitializeAntennaConstraint(*antenna_constraint,
                                    settings.antenna_constraint, antenna_names);
      }
    }

#ifdef HAVE_ARMADILLO
    ScreenConstraint* screen_constraint =
        dynamic_cast<ScreenConstraint*>(constraint.get());
    if (screen_constraint) {
      InitializeScreenConstraint(*screen_constraint,
                                 settings.screen_core_constraint,
                                 antenna_positions, source_positions);
    }
#endif

    SmoothnessConstraint* smoothness_constraint =
        dynamic_cast<SmoothnessConstraint*>(constraint.get());
    if (smoothness_constraint) {
      InitializeSmoothnessConstraint(*smoothness_constraint,
                                     settings.smoothness_ref_distance,
                                     antenna_positions);
    }
  }
}

}  // namespace ddecal
}  // namespace dp3