// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BDADDECAL_H_
#define DP3_BDADDECAL_H_

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

#include "../ddecal/Settings.h"
#include "../ddecal/SolutionWriter.h"
#include "../ddecal/gain_solvers/BdaSolverBuffer.h"
#include "../ddecal/gain_solvers/SolverBase.h"

#include "BDAResultStep.h"
#include "UVWFlagger.h"

namespace dp3 {
namespace steps {

/**
 * Direction-dependent calibration steps that supports Baseline Dependent
 * Averaging.
 *
 * BdaDdeCal internally has multiple substeps: For each direction, it creates
 * a series of steps. Each series starts with a Predict step. After the Predict
 * step, some optional steps may run. Finally, each series ends with a
 * BDAResultStep, which will contain the predicted visibilities for a direction.
 *
 * Workflow for BdaDdeCal:
 * 1. BdaDdeCal receives a BdaBuffer using its process() function. This buffer
 *    must contain visibilities and weights.
 *    It may also contain predicted model data from earlier BdaDdeCal steps.
 *    Those earlier steps should have applied the solutions to the model data.
 * 2. BdaDdeCal forwards a metadata-only copy of the BdaBuffer to the Predict
 *    steps for all directions.
 * 3. BdaDdeCal adds the predicted visiblities to its main input BdaBuffer.
 * 4.1. In only-predict mode, BdaDdeCal sums the predicted visibilities and
 *      sends them in a BdaBuffer to the next step. Processing is complete.
 * 4.2. Otherwise, BdaDdeCal appends its main input BdaBuffer to a
 *      BdaSolverBuffer.
 * 5. BdaSolverBuffer weights all input and model visibilities.
 * 6. When the BdaSolverBuffer has a complete solution interval, BdaDdeCal
 *    runs the solver on this interval. The solver uses the weighted data.
 * 7. If needed, the BdaSolverBuffer applies the solutions to the unweighted
 *    model data in the current solution interval.
 * 8. When the BdaSolverBuffer has a complete BdaBuffer:
 * 8.1. If subtracting is enabled, BdaDdeCal subtracts the model data (with
 *      the solutions applied) for all directions from the input data.
 * 8.2. If keeping model data is disabled, BdaDdeCal removes the model data
 *      for all directions from the BdaBuffer.
 * 8.3. BdaDdeCal passes the BdaBuffer to its next step. If keeping model data
 *      is enabled, this BdaBuffer includes the unweighted model data with
 *      the solutions applied.
 *
 * If both subtracting and keeping model data are disabled, the unweighted
 * model data is not needed anymore after weighting it. The BdaSolverBuffer
 * therefore only stores the weighted model data in that case, saving memory.
 *
 * The BdaBuffer structure/metadata remains equal for all BdaBuffers that are
 * involved in one iteration.
 *
 * The BdaSolverBuffer provides the functionality for solution intervals:
 * - BdaSolverBuffer detects when a solution interval is complete.
 * - BdaSolverBuffer only exposes the BDA rows in the interval to the solver.
 * - BdaSolverBuffer detects when a BdaBuffer is fully processed.
 */
class BdaDdeCal : public Step {
 public:
  /**
   * Constructor.
   * @param parset A parameter set with settings for this class.
   * @param prefix Prefix for reading settings from the parameter set.
   */
  BdaDdeCal(const common::ParameterSet& parset, const std::string& prefix);

  common::Fields getRequiredFields() const override;

  common::Fields getProvidedFields() const override {
    return (settings_.subtract ||
            (settings_.only_predict && !settings_.keep_model_data))
               ? kDataField
               : common::Fields();
  }

  bool process(std::unique_ptr<base::BdaBuffer>) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

  MsType outputs() const override { return MsType::kBda; };

  /**
   * Get the channel block index for a given channel.
   * @return Index of the channel block.
   */
  size_t GetChanBlockIndex(size_t channel, size_t n_channels,
                           size_t n_channel_blocks) const;

 private:
  void InitializeModelReuse();

  void InitializePredictSteps(const common::ParameterSet& parset,
                              const std::string& prefix);

  /// Initialize chan_block_start_freqs_.
  void DetermineChannelBlocks();

  /// Creates the list of directions for each sub-step and reused model data.
  void SetSourceDirections();

  /// @return A list with the center frequency for each channel block.
  std::vector<double> GetChannelBlockFrequencies() const;

  /// Extracts results from all sub-steps and appends them to model_buffers_.
  void ExtractResults();

  /// Checks if a buffer contains all named directions.
  bool HasAllDirections(const base::BdaBuffer& buffer) const;

  /// Sum all predictions directions into the main data buffer.
  void SumModels(base::BdaBuffer& buffer) const;

  /// Processes the data for the directions where all sub-steps gave results.
  void ProcessCompleteDirections();

  /// Solve the current solution interval using the BdaSolverBuffer
  /// and advance to the next solution interval.
  void SolveCurrentInterval();

  void InitializeCurrentSolutions();

  /// Extract the done buffers from the BdaSolverBuffer, subtract or delete
  /// model data if needed, and forward the buffers to the next step.
  void ProcessDone();

  void WriteSolutions();

 private:
  ddecal::Settings settings_;
  std::unique_ptr<ddecal::SolutionWriter> solution_writer_;

  /** The list of directions for the model data. */
  std::vector<base::Direction> source_directions_;

  /**
   * The first step of each sub-step chain.
   * There's a sub-step chain for each direction that does not reuse model data.
   */
  std::vector<std::shared_ptr<ModelDataStep>> steps_;
  /** Result step of each sub-step chain. */
  std::vector<std::shared_ptr<BDAResultStep>> result_steps_;

  std::unique_ptr<UVWFlagger> uvw_flagger_step_;
  std::shared_ptr<BDAResultStep> uvw_flagger_result_step_;

  /** For each direction, a list of patch names. */
  std::vector<std::vector<std::string>> patches_;

  /**
   * For each direction, the name for the model data in BdaBuffer.
   * First, there's an element for each sub-step chain.
   * Then, there's an element for each reused model data buffer.
   */
  std::vector<std::string> direction_names_;

  /** Stores the data buffers received from the process() function. */
  std::deque<std::unique_ptr<base::BdaBuffer>> input_buffers_;

  /**
   * Solver buffer.
   * Each process() call adds BdaBuffers to the solver buffer.
   * When the solution interval is complete, BdaDdeCal runs the solver and
   * removes old BdaBuffers from the solver buffer.
   * This variable is not used when only_predict is true.
   */
  std::unique_ptr<ddecal::BdaSolverBuffer> solver_buffer_;

  std::unique_ptr<ddecal::SolverBase> solver_;

  /** The solution interval, in seconds. */
  double solution_interval_duration_;

  /**
   * For each channel block, the start and end frequencies. The start and end
   * frequencies for channel block 'cb' are at indices 'cb' and 'cb + 1'.
   */
  std::vector<double> chan_block_start_freqs_;

  /**
   * Antenna lists, which contain the used antenna index for each baseline.
   * @{
   */
  std::vector<int> antennas1_;
  std::vector<int> antennas2_;
  /** @} */

  /// For each time, for each channel block, a vector of size nAntennas *
  /// nDirections
  std::vector<std::vector<std::vector<casacore::DComplex>>> solutions_;

  /// For each solution interval, the amount the solver iterations.
  std::vector<size_t> iterations_;
  std::vector<size_t> approx_iterations_;

  /// For each time, for each constraint, a vector of results (e.g. tec and
  /// phase)
  std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>
      constraint_solutions_;

  common::NSTimer timer_;
  common::NSTimer predict_timer_;
  common::NSTimer solve_timer_;
  common::NSTimer write_timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
