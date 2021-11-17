// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BDADDECAL_H
#define DP3_BDADDECAL_H

#include "Step.h"
#include "../ddecal/Settings.h"
#include "../ddecal/SolutionWriter.h"
#include "../ddecal/gain_solvers/BdaSolverBuffer.h"
#include "../ddecal/gain_solvers/SolverBase.h"

namespace dp3 {
namespace steps {

class BDAResultStep;
class ModelDataStep;

/**
 * Direction-dependent calibration steps that supports Baseline Dependent
 * Averaging.
 *
 * BdaDdeCal internally has multiple substeps: For each direction, it creates
 * a series of steps. Each series starts with a Predict step. After the Predict
 * step, some optional steps may run. Finally, each series ends with a
 * BDAResultStep, which will contain the predicted visibilities for a direction.
 *
 * Workflow for each iteration:
 * 1. BdaDdeCal receives a BDABuffer using its process() function. This buffer
 *    must contain visibilities and weights.
 * 2. BdaDdeCal forwards a metadata-only copy of the BDABuffer to the Predict
 *    steps for all directions.
 * 3.1. In only-predict mode, BdaDdeCal sums the predicted visibilities and
 *      sends them in a BDABuffer to the next step. Processing is complete.
 * 3.2. Otherwise, BdaDdeCal appends the received visibilities and weights and
 *      the predictions from the BDAResultSteps into a BdaSolverBuffer.
 * 4. When the BdaSolverBuffer has a complete solution interval, BdaDdeCal
 *    runs the solver on this interval.
 * 5. When subtracting is enabled, the solver applies the solutions to the
 *    predicted visibilities for each direction.
 * 6. When the BdaSolverBuffer has a complete BDABuffer:
 * 6.1. If subtracting is enabled, BdaDdeCal subtracts the predicted
 *      visibilities for all directions from the input data.
 * 6.2. BdaDdeCal passes the input data buffer to its next step.
 *
 * The BDABuffer structure/metadata remains equal for all BDABuffers that are
 * involved in one iteration.
 *
 * The BdaSolverBuffer provides the functionality for solution intervals:
 * - BdaSolverBuffer detects when a solution interval is complete.
 * - BdaSolverBuffer only exposes the BDA rows in the interval to the solver.
 * - BdaSolverBuffer detects when a BDABuffer is fully processed.
 */
class BdaDdeCal : public Step {
 public:
  /**
   * Constructor.
   * @param input_step Input step, for reading extra data.
   * @param parset A parameter set with settings for this class.
   * @param prefix Prefix for reading settings from the parameter set.
   */
  BdaDdeCal(InputStep* input_step, const common::ParameterSet& parset,
            const std::string& prefix);

  bool process(std::unique_ptr<base::BDABuffer>) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

  MsType outputs() const override { return MsType::kBda; };

 private:
  void InitializePredictSteps(InputStep* input,
                              const common::ParameterSet& parset,
                              const string& prefix);

  /// Initialize chan_block_start_freqs_.
  void DetermineChannelBlocks();

  /// @return A list with the first direction of each sub-step.
  std::vector<base::Direction> GetSourceDirections() const;

  /// @return A list with the center frequency for each channel block.
  std::vector<double> GetChannelBlockFrequencies() const;

  /// Extracts results from all sub-steps and appends them to model_buffers_.
  void ExtractResults();

  /// Processes the data for the directions where all sub-steps gave results.
  void ProcessCompleteDirections();

  /// Solve the current solution interval using the BdaSolverBuffer
  /// and advance to the next solution interval.
  void SolveCurrentInterval();

  void InitializeCurrentSolutions();

  void WriteSolutions();

 private:
  ddecal::Settings settings_;
  std::unique_ptr<ddecal::SolutionWriter> solution_writer_;

  /// For each direction, the first step of that direction.
  std::vector<std::shared_ptr<ModelDataStep>> steps_;
  /// For each direction, a result step.
  std::vector<std::shared_ptr<BDAResultStep>> result_steps_;
  /// For each direction, a list of patch names.
  std::vector<std::vector<std::string>> patches_;

  /**
   * Stores the data buffers received from the process() function.
   */
  std::deque<std::unique_ptr<base::BDABuffer>> input_buffers_;

  /**
   * Stores the result buffers from the sub steps. Each queue item holds
   * a vector of results for the directions.
   * When all directions ina queue item are valid, further processing can start.
   */
  std::deque<std::vector<std::unique_ptr<base::BDABuffer>>> model_buffers_;

  /**
   * Solver buffer.
   * Each process() call adds BDABuffers to the solver buffer.
   * When the solution interval is complete, BdaDdeCal runs the solver and
   * removes old BDABuffers from the solver buffer.
   * This variable is not used when only_predict is true.
   */
  std::unique_ptr<ddecal::BdaSolverBuffer> solver_buffer_;

  std::unique_ptr<ddecal::SolverBase> solver_;

  /** The solution interval, in seconds. */
  double solution_interval_;

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