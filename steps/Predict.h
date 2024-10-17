// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_PREDICT_H_
#define DP3_STEPS_PREDICT_H_

#include <dp3/steps/Step.h>

#include <mutex>

namespace dp3 {
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
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param input_type Input type, Regular (default) or Bda.
   */
  Predict(const common::ParameterSet& parset, const std::string& prefix,
          MsType input_type = MsType::kRegular);

  /**
   * Constructs the object with explicit source patterns.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param source_patterns Source patterns.
   * @param input_type Input type, Regular (default) or Bda.
   */
  Predict(const common::ParameterSet& parset, const std::string& prefix,
          const std::vector<std::string>& source_patterns,
          MsType input_type = MsType::kRegular);

  ~Predict() override {}

  common::Fields getRequiredFields() const override {
    // The actual work occurs in sub-steps, which come after Predict in the
    // step list.
    return {};
  }

  common::Fields getProvidedFields() const override {
    // The actual work occurs in sub-steps, which come after Predict in the
    // step list.
    return {};
  }

  bool process(std::unique_ptr<base::DPBuffer>) override;

  bool process(std::unique_ptr<base::BDABuffer>) override;

  void finish() override;

  void show(std::ostream&) const override;

  /// Ensures that all steps, including internal sub-steps, form a single list.
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
  void SetThreadData(std::mutex* mutex);

 private:
  /**
   * Common part of the constructors.
   * Parses parset arguments and sets up first_step_ and last_step_.
   */
  void Initialize(const common::ParameterSet& parset, const std::string& prefix,
                  MsType input_type);

  /// Input and output measurement set type: Regular (default) or Bda
  const MsType ms_type_;

  std::vector<std::shared_ptr<Step>> internal_steps_;
  std::shared_ptr<BDAAverager> bda_averager_;
  std::shared_ptr<OnePredict> predict_step_;
};

}  // namespace steps
}  // namespace dp3

#endif
