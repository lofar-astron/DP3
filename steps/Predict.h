// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_PREDICT_H
#define DP3_PREDICT_H

#include "Step.h"

namespace dp3 {
namespace base {
class PredictBuffer;
}
namespace common {
class ParameterSet;
}

namespace steps {
class OnePredict;

/**
 * @brief DP3 step class that predicts visibilities from a source model.
 * This step contains OnePredict sub-steps, that do the actual predictions and
 * optional pre- and postprocessing sub-steps for each OnePredict sub-step.
 */
class Predict : public ModelDataStep {
 public:
  /**
   * Constructs the object.
   * @param input_step Input step, for reading extra data.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   */
  Predict(InputStep& input_step, const common::ParameterSet& parset,
          const std::string& prefix);

  /**
   * Constructs the object with explicit source patterns.
   * @param input_step Input step, for reading extra data.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param source_patterns Source patterns.
   */
  Predict(InputStep& input_step, const common::ParameterSet& parset,
          const std::string& prefix,
          const std::vector<std::string>& source_patterns);

  virtual ~Predict() {}

  bool process(const base::DPBuffer& buffer) override;

  void finish() override;

  void show(std::ostream&) const override;

  void setNextStep(std::shared_ptr<Step> next_step) override;

  std::pair<double, double> GetFirstDirection() const override;

  void SetOperation(const std::string& operation);

  /**
   * Forwards thread synchronization structures to its predict sub-step.
   * @see OnePredict::SetThreadData().
   */
  void SetThreadData(aocommon::ThreadPool& pool, std::mutex* mutex);

  void SetPredictBuffer(std::shared_ptr<base::PredictBuffer> predict_buffer);

 private:
  std::shared_ptr<OnePredict> predict_step_;
};

}  // namespace steps
}  // namespace dp3

#endif
