// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverFactory.h"

#include <iostream>

#include "Settings.h"

#include "gain_solvers/DiagonalLowRankSolver.h"
#include "gain_solvers/DiagonalSolver.h"
#include "gain_solvers/FullJonesSolver.h"
#include "gain_solvers/LBFGSSolver.h"
#include "gain_solvers/HybridSolver.h"
#include "gain_solvers/IterativeDiagonalSolver.h"
#if defined(HAVE_CUDA_SOLVER)
#include "gain_solvers/IterativeDiagonalSolverCuda.h"
#endif
#include "gain_solvers/IterativeFullJonesSolver.h"
#include "gain_solvers/IterativeScalarSolver.h"
#include "gain_solvers/ScalarSolver.h"

#include "constraints/AmplitudeOnlyConstraint.h"
#include "constraints/AntennaConstraint.h"
#include "constraints/RotationConstraint.h"
#include "constraints/RotationAndDiagonalConstraint.h"
#ifdef ENABLE_SCREENFITTER
#include "constraints/ScreenConstraint.h"
#endif
#include "constraints/SmoothnessConstraint.h"
#include "constraints/TECConstraint.h"

namespace dp3 {
namespace ddecal {

namespace {

std::unique_ptr<SolverBase> CreateScalarSolver(SolverAlgorithm algorithm,
                                               const Settings& settings) {
  switch (algorithm) {
    case SolverAlgorithm::kDirectionIterative:
      return std::make_unique<IterativeScalarSolver>();
    case SolverAlgorithm::kDirectionSolve:
      return std::make_unique<ScalarSolver>();
    case SolverAlgorithm::kLowRank:
      throw std::runtime_error(
          "Low-rank algorithm only supports diagonal solving");
    case SolverAlgorithm::kLBFGS:
#ifdef HAVE_LIBDIRAC
      return std::make_unique<LBFGSSolver>(
          settings.lbfgs_robust_nu, settings.lbfgs_max_iter,
          settings.lbfgs_history_size, settings.lbfgs_minibatches,
          LBFGSSolver::kScalar);
#else
      throw std::runtime_error("LIBDIRAC for LBFGS not built");
#endif
    case SolverAlgorithm::kHybrid:
      break;  // CreateSolver should have handled this case.
  }
  assert(false);
  return nullptr;
}

std::unique_ptr<SolverBase> CreateDiagonalSolver(SolverAlgorithm algorithm,
                                                 const Settings& settings) {
#if defined(HAVE_CUDA_SOLVER)
  if (settings.use_gpu) {
    switch (algorithm) {
      case SolverAlgorithm::kDirectionIterative:
        return std::make_unique<IterativeDiagonalSolverCuda>(
            settings.keep_host_buffers);
      default:
        throw std::runtime_error(
            "usegpu=true, but no GPU implementation for solver algorithm is "
            "available.");
    }
  }
#else
  if (settings.use_gpu) {
    throw std::runtime_error(
        "usegpu=true, but DP3 is built without CUDA support.");
  }
#endif
  switch (algorithm) {
    case SolverAlgorithm::kDirectionIterative:
      return std::make_unique<IterativeDiagonalSolver>();
    case SolverAlgorithm::kDirectionSolve:
      return std::make_unique<DiagonalSolver>();
    case SolverAlgorithm::kLowRank:
      return std::make_unique<DiagonalLowRankSolver>(
          settings.n_lra_iterations, settings.n_lra_power_iterations);
    case SolverAlgorithm::kLBFGS:
#ifdef HAVE_LIBDIRAC
      return std::make_unique<LBFGSSolver>(
          settings.lbfgs_robust_nu, settings.lbfgs_max_iter,
          settings.lbfgs_history_size, settings.lbfgs_minibatches,
          LBFGSSolver::kDiagonal);
#else
      throw std::runtime_error("LIBDIRAC for LBFGS not built");
#endif
    case SolverAlgorithm::kHybrid:
      break;  // CreateSolver should have handled this case
  }
  assert(false);
  return nullptr;
}

std::unique_ptr<SolverBase> CreateFullJonesSolver(SolverAlgorithm algorithm,
                                                  const Settings& settings) {
  switch (algorithm) {
    case SolverAlgorithm::kDirectionIterative:
      return std::make_unique<IterativeFullJonesSolver>();
    case SolverAlgorithm::kDirectionSolve:
      return std::make_unique<FullJonesSolver>();
    case SolverAlgorithm::kLowRank:
      throw std::runtime_error(
          "Low-rank algorithm only supports diagonal solving");
    case SolverAlgorithm::kLBFGS:
#ifdef HAVE_LIBDIRAC
      return std::make_unique<LBFGSSolver>(
          settings.lbfgs_robust_nu, settings.lbfgs_max_iter,
          settings.lbfgs_history_size, settings.lbfgs_minibatches,
          LBFGSSolver::kFull);
#else
      throw std::runtime_error("LIBDIRAC for LBFGS not built");
#endif
    case SolverAlgorithm::kHybrid:
      break;  // CreateSolver should have handled this case.
  }
  assert(false);
  return nullptr;
}

void AddConstraints(SolverBase& solver, const Settings& settings,
                    const common::ParameterSet& parset,
                    const std::string& prefix) {
  if (settings.core_constraint != 0.0 || !settings.antenna_constraint.empty()) {
    solver.AddConstraint(std::make_unique<AntennaConstraint>());
  }
  if (settings.smoothness_constraint != 0.0) {
    solver.AddConstraint(std::make_unique<SmoothnessConstraint>(
        settings.smoothness_constraint, settings.smoothness_ref_frequency));
  }

  switch (settings.mode) {
    case base::CalType::kScalar:
    case base::CalType::kDiagonal:
    case base::CalType::kFullJones:
    case base::CalType::kScalarPhase:
    case base::CalType::kDiagonalPhase:
      // no extra constraints
      break;
    case base::CalType::kScalarAmplitude:
    case base::CalType::kDiagonalAmplitude:
      solver.AddConstraint(std::make_unique<AmplitudeOnlyConstraint>());
      break;
    case base::CalType::kTec:
    case base::CalType::kTecAndPhase: {
      const auto tec_mode = (settings.mode == base::CalType::kTec)
                                ? TECConstraint::Mode::kTecOnly
                                : TECConstraint::Mode::kTecAndCommonScalar;
      std::unique_ptr<TECConstraint> constraint;

      if (settings.approximate_tec) {
        auto approxConstraint =
            std::make_unique<ApproximateTECConstraint>(tec_mode);
        approxConstraint->SetMaxApproximatingIterations(
            settings.max_approx_iterations);
        approxConstraint->SetFittingChunkSize(settings.approx_chunk_size);
        constraint = std::move(approxConstraint);
      } else {
        constraint = std::make_unique<TECConstraint>(tec_mode);
      }
      constraint->setDoPhaseReference(settings.phase_reference);
      solver.AddConstraint(std::move(constraint));
      break;
    }
#ifdef ENABLE_SCREENFITTER
    case base::CalType::kTecScreen:
      solver.AddConstraint(
          std::make_unique<ScreenConstraint>(parset, prefix + "tecscreen."));
      break;
#endif
    case base::CalType::kRotationAndDiagonal: {
      auto constraint = std::make_unique<RotationAndDiagonalConstraint>();
      constraint->SetDoRotationReference(settings.rotation_reference);
      solver.AddConstraint(std::move(constraint));
      break;
    }
    case base::CalType::kRotation:
      solver.AddConstraint(std::make_unique<RotationConstraint>());
      break;
    default:
      throw std::runtime_error("Unexpected solving mode: " +
                               ToString(settings.mode));
  }
}

void InitializeSolver(SolverBase& solver, const Settings& settings) {
  solver.SetLLSSolverType(settings.lls_solver_type);
  solver.SetMaxIterations(settings.max_iterations);
  solver.SetAccuracy(settings.tolerance);
  solver.SetConstraintAccuracy(settings.approx_tolerance);
  solver.SetStepSize(settings.step_size);
  solver.SetDetectStalling(settings.detect_stalling, settings.step_diff_sigma);
}

std::unique_ptr<SolverBase> CreateSolver(const Settings& settings,
                                         const common::ParameterSet& parset,
                                         const std::string& prefix,
                                         SolverAlgorithm algorithm) {
  std::unique_ptr<SolverBase> solver;
  switch (settings.mode) {
    case base::CalType::kScalar:
    case base::CalType::kScalarAmplitude:
      solver = CreateScalarSolver(algorithm, settings);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kScalarPhase:
    case base::CalType::kTec:
    case base::CalType::kTecAndPhase:
      solver = CreateScalarSolver(algorithm, settings);
      solver->SetPhaseOnly(true);
      break;
    case base::CalType::kDiagonal:
    case base::CalType::kDiagonalAmplitude:
      solver = CreateDiagonalSolver(algorithm, settings);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kDiagonalPhase:
      solver = CreateDiagonalSolver(algorithm, settings);
      solver->SetPhaseOnly(true);
      break;
    case base::CalType::kFullJones:
    case base::CalType::kRotationAndDiagonal:
    case base::CalType::kRotation:
      solver = CreateFullJonesSolver(algorithm, settings);
      solver->SetPhaseOnly(false);
      break;
    case base::CalType::kTecScreen:
#ifdef ENABLE_SCREENFITTER
      solver = CreateScalarSolver(algorithm, settings);
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
        std::cout << "Warning: antenna constraint specifies antenna "
                  << constraint_name
                  << " which is not in the MS, ignoring it\n";
      } else {
        constraint_list.back().insert(constraint_list.back().end(),
                                      iter - antenna_names.begin());
      }
    }
  }
  constraint.SetAntennaSets(std::move(constraint_list));
}

#ifdef ENABLE_SCREENFITTER
void InitializeScreenConstraint(
    ScreenConstraint& constraint, double core_constraint,
    const std::vector<std::array<double, 3>>& antenna_positions,
    const std::vector<base::Direction>& source_directions) {
  constraint.InitPiercePoints(antenna_positions, source_directions);
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

std::unique_ptr<SolverBase> CreateSolver(const Settings& settings,
                                         const common::ParameterSet& parset,
                                         const std::string& prefix) {
  std::unique_ptr<SolverBase> solver;

  if (settings.solver_algorithm == SolverAlgorithm::kHybrid) {
    std::unique_ptr<SolverBase> a = CreateSolver(
        settings, parset, prefix, SolverAlgorithm::kDirectionSolve);
    // The max_iterations is divided by 6 to use at most 1/6th of the iterations
    // in the first solver.
    a->SetMaxIterations(std::max<size_t>(1u, settings.max_iterations / 6u));
    std::unique_ptr<SolverBase> b = CreateSolver(
        settings, parset, prefix, SolverAlgorithm::kDirectionIterative);

    auto hybrid_solver = std::make_unique<HybridSolver>();
    hybrid_solver->SetMaxIterations(settings.max_iterations);
    hybrid_solver->AddSolver(std::move(a));
    hybrid_solver->AddSolver(std::move(b));
    solver = std::move(hybrid_solver);
  } else {
    solver = CreateSolver(settings, parset, prefix, settings.solver_algorithm);
  }
  return solver;
}

void InitializeSolverConstraints(
    SolverBase& solver, const Settings& settings,
    const std::vector<std::array<double, 3>>& antenna_positions,
    const std::vector<std::string>& antenna_names,
    const std::vector<size_t>& solutions_per_direction,
    const std::vector<base::Direction>& source_directions,
    const std::vector<double>& frequencies) {
  for (const std::unique_ptr<Constraint>& constraint :
       solver.GetConstraints()) {
    // Initialize the constraint with some common metadata.
    constraint->Initialize(
        antenna_positions.size(),
        std::vector<uint32_t>(solutions_per_direction.begin(),
                              solutions_per_direction.end()),
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

#ifdef ENABLE_SCREENFITTER
    ScreenConstraint* screen_constraint =
        dynamic_cast<ScreenConstraint*>(constraint.get());
    if (screen_constraint) {
      InitializeScreenConstraint(*screen_constraint,
                                 settings.screen_core_constraint,
                                 antenna_positions, source_directions);
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
