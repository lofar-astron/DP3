// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Settings.h"
#include "../base/CalType.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <numeric>
#include <sstream>

#include <boost/algorithm/string/case_conv.hpp>

#include <aocommon/logger.h>

using dp3::base::CalType;

namespace dp3 {
namespace ddecal {

namespace {
std::string CreateParsetString(const dp3::common::ParameterSet& parset) {
  std::stringstream ss;
  ss << parset;
  return ss.str();
}

SolverAlgorithm ParseSolverAlgorithm(const std::string& str) {
  const std::string lowercase = boost::to_lower_copy(str);
  if (lowercase == "lowrank")
    return SolverAlgorithm::kLowRank;
  else if (lowercase == "directionsolve")
    return SolverAlgorithm::kDirectionSolve;
  else if (lowercase == "directioniterative")
    return SolverAlgorithm::kDirectionIterative;
  else if (lowercase == "hybrid")
    return SolverAlgorithm::kHybrid;
  else if (lowercase == "lbfgs")
    return SolverAlgorithm::kLBFGS;
  else
    throw std::runtime_error("Unknown solver algorithm specified: " + str);
}

SolverDataUse ParseSolverDataUse(const std::string& data_use_string) {
  const std::string lowercase = boost::to_lower_copy(data_use_string);
  if (lowercase == "single")
    return SolverDataUse::kSingle;
  else if (lowercase == "dual")
    return SolverDataUse::kDual;
  else if (lowercase == "full")
    return SolverDataUse::kFull;
  else
    throw std::runtime_error(
        "Unknown value specified for ddecal's parameter 'datause': " +
        data_use_string + ". Please use one of single, dual, full");
}

}  // namespace

std::string ToString(SolverAlgorithm algorithm) {
  switch (algorithm) {
    case SolverAlgorithm::kLowRank:
      return "lowrank";
    case SolverAlgorithm::kDirectionSolve:
      return "directionsolve";
    case SolverAlgorithm::kDirectionIterative:
      return "directioniterative";
    case SolverAlgorithm::kHybrid:
      return "hybrid";
    case SolverAlgorithm::kLBFGS:
      return "LBFGS";
  }
  return "invalid algorithm";
}

Settings::Settings(const common::ParameterSet& _parset,
                   const std::string& _prefix)
    : parset(&_parset),
      name(_prefix),
      h5parm_name(parset->isDefined(_prefix + "h5parm")
                      ? GetString("h5parm")
                      : parset->getString("msin") + "/instrument.h5"),
      stat_filename(GetString("statfilename", "")),
      parset_string(CreateParsetString(_parset)),
      mode(dp3::base::StringToCalType(
          boost::to_lower_copy(GetString("mode", "complexgain")))),
      propagate_solutions(GetBool("propagatesolutions", false)),
      propagate_converged_only(GetBool("propagateconvergedonly", false)),
      flag_unconverged(GetBool("flagunconverged", false)),
      flag_diverged_only(GetBool("flagdivergedonly", false)),
      only_predict(GetBool("onlypredict", false)),
      subtract(GetBool("subtract", false)),
      keep_model_data(GetBool("keepmodel", false)),
      solver_algorithm(
          ParseSolverAlgorithm(GetString("solveralgorithm", "directionsolve"))),
      solution_interval(GetUint("solint", 1)),
      min_vis_ratio(GetDouble("minvisratio", 0.0)),
      n_channels(GetUint("nchan", 1)),
      solutions_per_direction(GetSizeTVector("solutions_per_direction", {})),
      // Constraints
      model_weighted_constraints(GetBool("model_weighted_constraints", false)),
      core_constraint(GetDouble("coreconstraint", 0.0)),
      antenna_constraint(ReadAntennaConstraint()),
      smoothness_constraint(GetDouble("smoothnessconstraint", 0.0)),
      smoothness_ref_frequency(GetDouble("smoothnessreffrequency", 0.0)),
      smoothness_ref_distance(GetDouble("smoothnessrefdistance", 0.0)),
      smoothness_spectral_exponent(
          GetDouble("smoothnessspectralexponent", -1.0)),
      smoothness_kernel_truncation(
          GetBool("smoothness_kernel_truncation", true)),
      smoothness_dd_factors(GetDoubleVector("smoothness_dd_factors")),
      screen_core_constraint(GetDouble("tecscreen.coreconstraint", 0.0)),

      // Solver settings
      lls_solver_type(
          ddecal::LLSSolver::ParseType(GetString("llssolver", "qr"))),
      max_iterations(GetUint("maxiter", 50)),
      tolerance(GetDouble("tolerance", 1.0e-5)),
      step_size(GetDouble("stepsize", 0.2)),
      solver_data_use(ParseSolverDataUse(GetString("datause", "full"))),
      detect_stalling(GetBool("detectstalling", true)),
      step_diff_sigma(detect_stalling ? GetDouble("stepsigma", 0.1) : 0.1),
      // Only read these settings when needed: If it is defined, but not used,
      // the application will give a warning.
      approximate_tec((mode == CalType::kTec || mode == CalType::kTecAndPhase)
                          ? GetBool("approximatetec", false)
                          : false),
      phase_reference((mode == CalType::kTec || mode == CalType::kTecAndPhase)
                          ? GetBool("phasereference", true)
                          : false),
      approx_tolerance(approximate_tec
                           ? GetDouble("approxtolerance", tolerance * 10.0)
                           : 0.0),
      max_approx_iterations(
          approximate_tec ? GetUint("maxapproxiter", max_iterations / 2) : 0),
      approx_chunk_size(approximate_tec ? GetUint("approxchunksize", 0) : 0),
      rotation_reference((mode == CalType::kRotationAndDiagonal)
                             ? GetBool("rotationreference", false)
                             : false),
      rotation_diagonal_mode(
          (mode == CalType::kRotationAndDiagonal)
              ? dp3::base::StringToCalType(boost::to_lower_copy(
                    GetString("rotationdiagonalmode", "diagonal")))
              : CalType::kDiagonal),
      lbfgs_robust_nu((solver_algorithm == SolverAlgorithm::kLBFGS)
                          ? GetDouble("solverlbfgs.dof", 200.0)
                          : 200.0),
      lbfgs_max_iter((solver_algorithm == SolverAlgorithm::kLBFGS)
                         ? GetUint("solverlbfgs.iter", 4)
                         : 4),
      lbfgs_history_size((solver_algorithm == SolverAlgorithm::kLBFGS)
                             ? GetUint("solverlbfgs.history", 10)
                             : 10),
      lbfgs_minibatches((solver_algorithm == SolverAlgorithm::kLBFGS)
                            ? GetUint("solverlbfgs.minibatches", 1)
                            : 1),
      lbfgs_min_solution((solver_algorithm == SolverAlgorithm::kLBFGS)
                             ? GetDouble("solverlbfgs.min_solution", 0.0)
                             : 0.0),
      lbfgs_max_solution((solver_algorithm == SolverAlgorithm::kLBFGS)
                             ? GetDouble("solverlbfgs.max_solution", 0.0)
                             : 0.0),
      use_gpu(GetBool("usegpu", 0)),
      keep_host_buffers(GetBool("keep_host_buffers", 0)),
      n_lra_iterations((solver_algorithm == SolverAlgorithm::kLowRank)
                           ? GetUint("lra.iterations", 25)
                           : 1),
      n_lra_power_iterations((solver_algorithm == SolverAlgorithm::kLowRank)
                                 ? GetUint("lra.power_iterations", 10)
                                 : 1),
      // Column reader settings
      model_data_columns(ReadModelDataColumns()),
      reuse_model_data(GetStringVector("reusemodel")),

      // IDG settings
      idg_region_filename(GetString("idg.regions", "")),
      idg_image_filenames(GetStringVector("idg.images")),

      // Sagecal predict
      use_sagecal_predict(GetBool("sagecalpredict", false)),

      directions(GetStringVector("directions")),
      source_db(GetString("sourcedb", ""))

{
  for (double factor : smoothness_dd_factors) {
    // Factor is disallowed to be >1 because the size of the kernel is currently
    // not enlarged.
    if (factor <= 0.0 || factor > 1.0)
      throw std::runtime_error(
          "The values of smoothness_dd_factors should be larger than zero and "
          "at most one (i.e. they can only shrink the smoothing kernel)");
  }
}

void Settings::PrepareSolutionsPerDirection(size_t n_directions) {
  if (solutions_per_direction.size() > n_directions) {
    throw std::runtime_error(
        "The size of solutions_per_direction should be less or equal "
        "than the number of directions.");
  }

  // Pad itsSolutionsPerDirection with 1s
  solutions_per_direction.resize(n_directions, 1);

  if (std::find(solutions_per_direction.begin(), solutions_per_direction.end(),
                0) != solutions_per_direction.end()) {
    throw std::runtime_error(
        "All entries in solutions_per_direction should be > 0.");
  }

  for (size_t val : solutions_per_direction) {
    if (solution_interval % val != 0) {
      throw std::runtime_error(
          "Values in solutions_per_direction should be integer divisors "
          "of solint (" +
          std::to_string(solution_interval) + "), " + std::to_string(val) +
          " is not.");
    }
  }

  const size_t max_n_solutions_per_direction = *std::max_element(
      solutions_per_direction.begin(), solutions_per_direction.end());

  if (max_n_solutions_per_direction > 1) {
    // Since info().ntime() might not be set at this stage, throw an error
    // in case itsRequestedSolInt equals 0, and DD intervals are used
    if (solution_interval == 0) {
      throw std::runtime_error(
          "Can't combine direction-dependent solution intervals with solint=0. "
          "Either set solint to a non-zero value, or set all "
          "solutions_per_direction entries to 1.");
    }

    const size_t actual_solution_interval =
        solution_interval / max_n_solutions_per_direction;
    if (actual_solution_interval == 0) {
      throw std::runtime_error(
          "Maximum value in solutions_per_direction of ddecal settings is "
          "larger than solint value.");
    }
  }

  if (!smoothness_dd_factors.empty() &&
      smoothness_dd_factors.size() != n_directions) {
    throw std::runtime_error(
        "Invalid number of values specified for the direction-dependent "
        "smoothness factors. This number should be equal to the number of "
        "directions solved for.");
  }
}

size_t Settings::GetNSolutions() const {
  return std::accumulate(solutions_per_direction.begin(),
                         solutions_per_direction.end(), 0u);
}

bool Settings::GetBool(const std::string& key, bool default_value) const {
  return parset->getBool(name + key, default_value);
}

unsigned int Settings::GetUint(const std::string& key,
                               unsigned int default_value) const {
  return parset->getUint(name + key, default_value);
}

std::vector<size_t> Settings::GetSizeTVector(
    const std::string& key, const std::vector<size_t>& default_value) const {
  std::vector<unsigned int> uint_vector = parset->getUintVector(
      name + key,
      std::vector<unsigned int>(default_value.begin(), default_value.end()));
  return std::vector<size_t>(uint_vector.begin(), uint_vector.end());
}

double Settings::GetDouble(const std::string& key, double default_value) const {
  return parset->getDouble(name + key, default_value);
}

std::string Settings::GetString(const std::string& key) const {
  return parset->getString(name + key);
}

std::string Settings::GetString(const std::string& key,
                                const std::string default_value) const {
  return parset->getString(name + key, default_value);
}

std::vector<std::string> Settings::GetStringVector(
    const std::string& key) const {
  return parset->getStringVector(name + key, std::vector<std::string>());
}

std::vector<double> Settings::GetDoubleVector(const std::string& key) const {
  return parset->getDoubleVector(name + key, std::vector<double>());
}

std::vector<std::set<std::string>> Settings::ReadAntennaConstraint() const {
  std::vector<std::set<std::string>> antenna_constraint;

  const std::vector<std::string> constraint_list =
      GetStringVector("antennaconstraint");
  for (const std::string& constraint_name : constraint_list) {
    dp3::common::ParameterValue constraint_param(constraint_name);
    std::vector<std::string> antenna_list = constraint_param.getStringVector();
    antenna_constraint.emplace_back(antenna_list.begin(), antenna_list.end());
    // Check the size using the created set, and not using the vector, since
    // the set removes duplicate antenna names.
    if (antenna_constraint.back().size() == 1)
      throw std::runtime_error(
          "Error: antennaconstraint given that should constrain a group of "
          "antennas with one antenna in it. This does not make sense (did "
          "you forget using two square brackets? [[ ant1, ant2 ]] )");
  }

  return antenna_constraint;
}

std::vector<std::string> Settings::ReadModelDataColumns() const {
  std::vector<std::string> columns = GetStringVector("modeldatacolumns");

  // The statement below allows DDECal to be backwards compatible, e.g.
  // DP3 msin.modelcolumn=MY_MODEL_DATA ddecal.usemodelcolumn=true
  // msin=tDDECal.MS msout=.
  if (columns.empty() && GetBool("usemodelcolumn", false)) {
    aocommon::Logger::Warn
        << "Warning: The input contains the deprecated " + name +
               "usemodelcolumn setting, possibly combined with the "
               "deprecated msin.modelcolumn setting. Please use " +
               name + "modeldatacolumns instead.\n";
    columns.push_back(parset->getString("msin.modelcolumn", "MODEL_DATA"));
  }

  return columns;
}

std::vector<double> Settings::GetExpandedSmoothnessDdFactors() const {
  const size_t n_directions = solutions_per_direction.size();
  if (n_directions == 0 || smoothness_dd_factors.empty()) {
    return smoothness_dd_factors;
  } else {
    std::vector<double> result;
    for (size_t d = 0; d != n_directions; ++d) {
      const double direction_factor = smoothness_dd_factors[d];
      for (size_t i = 0; i != solutions_per_direction[d]; ++i) {
        result.emplace_back(direction_factor);
      }
    }
    return result;
  }
}

void ShowConstraintSettings(std::ostream& output, const Settings& settings) {
  using dp3::common::operator<<;
  if (!settings.antenna_constraint.empty())
    output << "  antennaconstraint:   " << settings.antenna_constraint << '\n';
  if (settings.core_constraint != 0.0)
    output << "  coreconstraint:      " << settings.core_constraint << '\n';
  if (settings.smoothness_constraint != 0.0)
    output << "  smoothnessconstraint:" << settings.smoothness_constraint
           << '\n';
  if (settings.smoothness_ref_frequency != 0.0)
    output << "  smoothnessreffrequency:" << settings.smoothness_ref_frequency
           << '\n';
  if (settings.smoothness_ref_distance != 0.0)
    output << "  smoothnessrefdistance:" << settings.smoothness_ref_distance
           << '\n';
  if (settings.screen_core_constraint != 0.0)
    output << "  tecscreen.coreconstraint:" << settings.screen_core_constraint
           << '\n';
}

std::vector<size_t> GetSolutionToDirectionVector(
    const std::vector<uint32_t>& solutions_per_direction) {
  std::vector<size_t> result;
  result.reserve(solutions_per_direction.size());
  for (size_t d = 0; d != solutions_per_direction.size(); ++d) {
    for (size_t i = 0; i != solutions_per_direction[d]; ++i) {
      result.emplace_back(d);
    }
  }
  return result;
}

}  // namespace ddecal
}  // namespace dp3
