// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BDADDECAL_H
#define DP3_BDADDECAL_H

#include "Step.h"
#include "../ddecal/Settings.h"
#include "../ddecal/gain_solvers/BDASolverBuffer.h"

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
 *      the predictions from the BDAResultSteps into a BDASolverBuffer.
 * 4. When the BDASolverBuffer has a complete solution interval, BdaDdeCal
 *    runs the solver on this interval.
 * 5. When subtracting is enabled, the solver applies the solutions to the
 *    predicted visibilities for each direction.
 * 6. When the BDASolverBuffer has a complete BDABuffer:
 * 6.1. If subtracting is enabled, BdaDdeCal subtracts the predicted
 *      visibilities for all directions from the input data.
 * 6.2. BdaDdeCal passes the input data buffer to its next step.
 *
 * The BDABuffer structure/metadata remains equal for all BDABuffers that are
 * involved in one iteration.
 *
 * The BDASolverBuffer provides the functionality for solution intervals:
 * - BDASolverBuffer detects when a solution interval is complete.
 * - BDASolverBuffer only exposes the BDA rows in the interval to the solver.
 * - BDASolverBuffer detects when a BDABuffer is fully processed.
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

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

  MsType outputs() const override { return MsType::kBda; };

 private:
  void InitializePredictSteps(InputStep* input,
                              const common::ParameterSet& parset,
                              const string& prefix);

 private:
  ddecal::Settings settings_;

  /// For each direction, the first step of that direction.
  std::vector<std::shared_ptr<ModelDataStep>> steps_;
  /// For each direction, a result step.
  std::vector<std::shared_ptr<BDAResultStep>> result_steps_;

  /**
   * Stores the data buffers received from the process() function.
   * When a buffer is fully processed, it sends it to the next step.
   * This member is not used when only_predict is true.
   */
  // std::deque<std::unique_ptr<base::BDABuffer>> data_buffers;

  /**
   * Solver buffer.
   * Each process() call adds BDABuffers to the solver buffer.
   * When the solution interval is complete, BdaDdeCal runs the solver and
   * removes old BDABuffers from the solver buffer.
   * This variable is not used when only_predict is true.
   */
  ddecal::BDASolverBuffer solver_buffer_;
};

}  // namespace steps
}  // namespace dp3

#endif