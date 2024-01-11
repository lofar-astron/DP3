// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Settings.h"
#include "../base/CalType.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <sstream>

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
      core_constraint(GetDouble("coreconstraint", 0.0)),
      antenna_constraint(ReadAntennaConstraint()),
      smoothness_constraint(GetDouble("smoothnessconstraint", 0.0)),
      smoothness_ref_frequency(GetDouble("smoothnessreffrequency", 0.0)),
      smoothness_ref_distance(GetDouble("smoothnessrefdistance", 0.0)),
      screen_core_constraint(GetDouble("tecscreen.coreconstraint", 0.0)),

      // Solver settings
      lls_solver_type(
          ddecal::LLSSolver::ParseType(GetString("llssolver", "qr"))),
      max_iterations(GetUint("maxiter", 50)),
      tolerance(GetDouble("tolerance", 1.e-4)),
      step_size(GetDouble("stepsize", 0.2)),
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
      source_db(GetString("sourcedb", "")) {}

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
}  // namespace ddecal
}  // namespace dp3
