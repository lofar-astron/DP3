// DDECal.h: DP3 step class to calibrate direction dependent gains
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_DDECAL_H_
#define DP3_STEPS_DDECAL_H_

#include <fstream>
#include <string>
#include <vector>

#include <aocommon/recursivefor.h>

#include <schaapcommon/h5parm/h5cache.h>
#include <schaapcommon/h5parm/jonesparameters.h>

#include "common/ParameterSet.h"

#include "ddecal/Settings.h"
#include "ddecal/SolutionWriter.h"
#include "ddecal/constraints/Constraint.h"
#include "ddecal/gain_solvers/SolverBase.h"

#include "MultiResultStep.h"
#include "ResultStep.h"
#include "UVWFlagger.h"

namespace dp3 {
namespace steps {

/// @brief This class is a Step class to calibrate (direction dependent) gains.
class DDECal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DDECal(const common::ParameterSet& parameterSet, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  common::Fields getProvidedFields() const override {
    return (settings_.subtract ||
            (settings_.only_predict && !settings_.keep_model_data))
               ? kDataField
               : common::Fields();
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void checkMinimumVisibilities(size_t bufferIndex);

  void flagChannelBlock(size_t cbIndex, size_t bufferIndex);

  /// Call the actual solver (called once per solution interval)
  void doSolve();

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

 private:
  void initializeColumnReaders(const common::ParameterSet&,
                               const std::string& prefix);
  void initializeModelReuse();
  void initializeInitialSolutionsH5Parm(const common::ParameterSet& parset,
                                        const std::string& prefix);
  void initializeIDG(const common::ParameterSet& parset,
                     const std::string& prefix);
  void initializePredictSteps(const common::ParameterSet& parset,
                              const std::string& prefix);

  void setModelNextSteps(Step&, const std::string& direction,
                         const common::ParameterSet& parset,
                         const std::string& prefix) const;

  void doPrepare();

  /// Initializes solutions for a new solution interval.
  /// Based on progation settings, either copies the previous solution or
  /// writes default values to the new solution.
  /// @param buffer_index Index within the current solution interval set.
  void InitializeSolutions(size_t buffer_index);

  /// Write all solutions to an H5Parm file using itsSolutionWriter.
  void WriteSolutions();

  /// Sums all model data buffers into the main data buffer.
  /// Removes the model data buffers if 'keepmodel' is false.
  /// This function implements the behavior for the "onlypredict" setting.
  void SumModels(size_t buffer_index);

  /// Applies a single solution to all directions.
  /// (Helper function for CorrectAndSubtractModels.)
  void ApplySolution(
      base::DPBuffer& buffer, size_t baseline, size_t channel,
      const std::vector<std::complex<double>>& channel_block_solutions) const;

  /// Applies the solutions to the model data for all directions.
  /// If "keepmodel" is true, overwrites the model data with the corrected model
  /// data. If "keepmodel" is false, removes the model data buffers.
  /// If "subtract" is true, subtracts all corrected model data from the main
  /// input data buffer.
  void CorrectAndSubtractModels(size_t buffer_index);

  /// Read the Jones matrix for a single time step and a single direction from
  /// one or two solution tables.
  xt::xtensor<std::complex<float>, 3> ReadJonesMatrixFromH5Parm(
      const base::Direction& direction, double timestamp,
      schaapcommon::h5parm::GainType gain_type,
      schaapcommon::h5parm::SolTab* first_soltab,
      schaapcommon::h5parm::SolTab* second_soltab);

  ddecal::Settings settings_;

  /// The input data buffers for the current set of solution intervals.
  /// Maximum dimensions: itsSolIntCount x itsRequestedSolInt
  std::vector<std::vector<std::unique_ptr<base::DPBuffer>>> input_buffers_;
  /// Original flags of the input buffers for the current solution interval.
  /// This member is only used if itsUVWFlagger is active.
  /// Dimensions: ( solution_interval x step_within_interval x baseline x
  /// channel x correlation )
  xt::xtensor<bool, 5> original_flags_;

  /// The time of the current buffer (in case of solint, average time)
  double average_time_;

  /// For each time, for each channel block, a vector of size nAntennas *
  /// SolverBase::NSolutions() * nPolarizations, with nPolarizations changing
  /// fastest.
  std::vector<std::vector<std::vector<casacore::DComplex>>> solutions_;
  std::vector<size_t> n_iterations_;  // Number of iterations taken
  std::vector<size_t> n_approximating_iterations_;

  /// For each time, for each constraint, a vector of results (e.g. tec and
  /// phase)
  std::vector<std::vector<std::vector<ddecal::ConstraintResult>>>
      constraint_solutions_;

  std::unique_ptr<ddecal::SolutionWriter> solution_writer_;

  /// Number of timeslots to store per solution interval as requested
  /// by the user in the parset.
  size_t requested_solution_interval_;
  size_t n_solution_intervals_;  ///< Number of solution intervals to buffer
  /// Index of the first solution in the current solution interval set.
  size_t first_solution_index_;
  size_t n_channels_;
  /// For each channel block, the nr of unflagged vis and the total nr of vis.
  std::vector<std::pair<size_t, size_t>> visibilities_in_interval_;
  /// For each channel block, the index in the channels at which this channel
  /// block starts.
  std::vector<size_t> channel_block_start_;
  std::vector<double> channel_block_frequencies_;
  /// For each direction, a vector of patches.
  std::vector<std::vector<std::string>> patches_per_direction_;
  /// For each direction, the name for the model data in DPBuffer.
  std::vector<std::string> direction_names_;
  /// Expanded version of reusemodel patterns.
  std::vector<std::string> reused_direction_names_;
  /// Maps direction indices to the cluster central direction.
  std::vector<base::Direction> source_directions_;

  /// First antenna for each baseline. Contains used antennas only.
  std::vector<int> antennas1_;
  /// Second antenna for each baseline. Contains used antennas only.
  std::vector<int> antennas2_;
  std::vector<double> weights_per_antenna_;

  UVWFlagger uvw_flag_step_;
  /// Result step for data after UV-flagging
  std::shared_ptr<ResultStep> data_result_step_;
  /// For each direction, the first step in the chain that computes the model.
  /// When reusing model data, the step for that direction is empty/null.
  std::vector<std::shared_ptr<ModelDataStep>> steps_;
  /// For each direction, the required fields of the step chain.
  std::vector<common::Fields> required_fields_;
  /// For each directions, a multiresultstep with all times.
  /// When reusing model data, the result step for that direction is empty/null.
  std::vector<std::shared_ptr<MultiResultStep>> result_steps_;

  /// Store the solution for later steps of processing in DPBuffer. Note: only
  /// works for 1 direction.
  bool store_solution_in_buffer_;

  /// Stores the H5Parm file and loads all solutions into memory when the user
  /// requests the solver to use initial solutions.
  std::unique_ptr<schaapcommon::h5parm::H5Parm> initial_solutions_;
  std::string initial_solutions_h5_parm_name;
  std::vector<std::string> initial_solutions_table_;
  std::vector<schaapcommon::h5parm::SolTab> solution_tables_;
  bool initial_solutions_are_full_jones_;
  /// Specifies the InterpolationType, MissingAntennaBehavior, and GainType for
  /// extracting the Jones parameters from itsInitialSolutions.
  /// @{
  schaapcommon::h5parm::JonesParameters::InterpolationType interpolation_type_;
  schaapcommon::h5parm::JonesParameters::MissingAntennaBehavior
      missing_antenna_behavior_;
  std::vector<schaapcommon::h5parm::GainType> gain_types_;
  /// @}

  common::NSTimer timer_;
  common::NSTimer predict_timer_;
  common::NSTimer solve_timer_;
  common::NSTimer write_timer_;
  std::mutex measures_mutex_;
  std::unique_ptr<ddecal::SolverBase> solver_;
  std::unique_ptr<std::ofstream> statistics_stream_;
};

}  // namespace steps
}  // namespace dp3

#endif
