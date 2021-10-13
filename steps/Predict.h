// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_PREDICT_H
#define DP3_PREDICT_H

#include "Step.h"

#include <mutex>

namespace aocommon {
class ThreadPool;
}

namespace dp3 {
namespace base {
class PredictBuffer;
}
namespace common {
class ParameterSet;
}

namespace steps {
class BDAAverager;
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
   * @param input_type Input type, Regular (default) or Bda.
   */
  Predict(InputStep& input_step, const common::ParameterSet& parset,
          const std::string& prefix, MsType input_type = MsType::kRegular);

  /**
   * Constructs the object with explicit source patterns.
   * @param input_step Input step, for reading extra data.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param source_patterns Source patterns.
   * @param input_type Input type, Regular (default) or Bda.
   */
  Predict(InputStep& input_step, const common::ParameterSet& parset,
          const std::string& prefix,
          const std::vector<std::string>& source_patterns,
          MsType input_type = MsType::kRegular);

  virtual ~Predict() {}

  bool process(const base::DPBuffer& buffer) override;

  bool process(std::unique_ptr<base::BDABuffer>) override;

  void finish() override;

  void show(std::ostream&) const override;

  void setNextStep(std::shared_ptr<Step> next_step) override;

  base::Direction GetFirstDirection() const override;

  void SetOperation(const std::string& operation);

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType dt) const override { return dt == ms_type_; }

  MsType outputs() const override { return ms_type_; }

  /**
   * Forwards thread synchronization structures to its predict sub-step.
   * @see OnePredict::SetThreadData().
   */
  void SetThreadData(aocommon::ThreadPool& pool, std::mutex* mutex);

  void SetPredictBuffer(std::shared_ptr<base::PredictBuffer> predict_buffer);

 private:
  /**
   * Common part of the constructors.
   * Parses parset arguments and sets up first_step_ and last_step_.
   */
  void Initialize(InputStep& input_step, const common::ParameterSet& parset,
                  const string& prefix, MsType input_type);

  /// Input and output measurement set type: Regular (default) or Bda
  const MsType ms_type_;

  std::vector<std::shared_ptr<Step>> steps_before_predict_;
  std::shared_ptr<BDAAverager> bda_averager_;
  std::shared_ptr<OnePredict> predict_step_;
  std::vector<std::shared_ptr<Step>> steps_after_predict_;
};

}  // namespace steps
}  // namespace dp3

#endif
