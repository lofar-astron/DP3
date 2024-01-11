// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SETTINGS_H
#define DP3_DDECAL_SETTINGS_H

#include "linear_solvers/LLSSolver.h"

#include <set>
#include <string>
#include <vector>

namespace dp3 {
namespace base {
enum class CalType;
}
namespace common {
class ParameterSet;
}

namespace ddecal {

enum class SolverAlgorithm {
  kLowRank,
  kDirectionSolve,
  kDirectionIterative,
  kHybrid,
  kLBFGS
};

std::string ToString(SolverAlgorithm algorithm);

/// @brief This struct parses the DDECal parset settings and stores them.
struct Settings {
 public:
  /**
   * Construct the object by reading settings from a parameter set.
   * @param parset A parameter set with DDECal settings.
   * @param prefix The prefix for accessing the parameter set.
   */
  Settings(const common::ParameterSet& parset, const std::string& prefix);

 private:
  /**
   * Retrieve an optional boolean from the parset.
   */
  bool GetBool(const std::string& key, bool default_value) const;

  /**
   * Retrieve an optional unsigned integer from the parset.
   */
  unsigned int GetUint(const std::string& key,
                       unsigned int default_value) const;

  /**
   * @brief Retrieve unsigned integer vector from the parset. If not found,
   * the default vector is returned.
   */
  std::vector<size_t> GetSizeTVector(
      const std::string& key, const std::vector<size_t>& default_value) const;

  /**
   * Retrieve an optional double from the parset.
   */
  double GetDouble(const std::string& key, double default_value) const;

  /**
   * Retrieve a mandatory string from the parset.
   */
  std::string GetString(const std::string& key) const;

  /**
   * Retrieve an optional string from the parset.
   */
  std::string GetString(const std::string& key,
                        const std::string default_value) const;

  /**
   * Retrieve an optional string vector from the parset.
   * @return If found: The found strings. If not found: Empty vector.
   */
  std::vector<std::string> GetStringVector(const std::string& key) const;

  std::vector<std::set<std::string>> ReadAntennaConstraint() const;

  std::vector<std::string> ReadModelDataColumns() const;

 private:
  const common::ParameterSet* parset;  // Only valid during construction!

 public:
  const std::string name;
  const std::string h5parm_name;
  const std::string stat_filename;
  const std::string parset_string;

  const base::CalType mode;
  const bool propagate_solutions;
  const bool propagate_converged_only;
  const bool flag_unconverged;
  const bool flag_diverged_only;
  const bool only_predict;
  const bool subtract;
  const bool keep_model_data;
  const SolverAlgorithm solver_algorithm;

  const size_t solution_interval;
  const double min_vis_ratio;
  const size_t n_channels;

  // Constraint settings.
  const double core_constraint;
  const std::vector<std::set<std::string>> antenna_constraint;
  const double smoothness_constraint;
  const double smoothness_ref_frequency;
  const double smoothness_ref_distance;
  const double screen_core_constraint;

  // Solver settings.
  const ddecal::LLSSolverType lls_solver_type;
  const size_t max_iterations;
  const double tolerance;
  const double step_size;
  const bool detect_stalling;
  const double step_diff_sigma;
  const bool approximate_tec;
  const bool phase_reference;
  const double approx_tolerance;
  const size_t max_approx_iterations;
  const size_t approx_chunk_size;
  const bool rotation_reference;
  // LBFGS robust parameter (aka degrees of freedom)
  const double lbfgs_robust_nu;
  // LBFGS max iterations per mini-batch
  const size_t lbfgs_max_iter;
  // LBFGS history size
  const size_t lbfgs_history_size;
  // LBFGS minibatches
  const size_t lbfgs_minibatches;
  const bool use_gpu;
  // keep host buffers between solve iteration
  // for the GPU solver
  const bool keep_host_buffers;
  // Number of iterations for the low-rank approximation (LRA) method
  const size_t n_lra_iterations;
  // In each lra iteration, the number of power-method iterations to take
  const size_t n_lra_power_iterations;

  const std::vector<std::string> model_data_columns;
  const std::vector<std::string> reuse_model_data;

  const std::string idg_region_filename;
  const std::vector<std::string> idg_image_filenames;

  const bool use_sagecal_predict;

  const std::vector<std::string> directions;
  const std::string source_db;
};

/** Writes the relevant constraints of the @a settings to the @a output. */
void ShowConstraintSettings(std::ostream& output, const Settings& settings);

}  // namespace ddecal
}  // namespace dp3

#endif
