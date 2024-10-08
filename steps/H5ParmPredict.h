// H5ParmPredict.h: DP3 step class to H5ParmPredict visibilities from a source
// model
//
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to H5ParmPredict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_H5PARM_PREDICT_H_
#define DP3_STEPS_H5PARM_PREDICT_H_

#include <dp3/base/DP3.h>
#include <dp3/steps/Step.h>

#include "../common/Timer.h"

#include "Predict.h"
#include "ResultStep.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step class to predict visibilities using an H5Parm file
/// with a source model.
class H5ParmPredict : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  H5ParmPredict(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    // Combine the result of all sub steps.
    return base::GetChainRequiredFields(itsPredictSteps.front());
  }

  common::Fields getProvidedFields() const override {
    // Combine the result of all sub steps.
    common::Fields fields;
    for (std::shared_ptr<Step> step = itsPredictSteps.front(); step;
         step = step->getNextStep()) {
      fields |= step->getProvidedFields();
    }
    return fields;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

 private:
  std::string itsName;

  std::vector<std::shared_ptr<Predict>> itsPredictSteps;
  std::shared_ptr<ResultStep> itsResultStep;

  std::string itsH5ParmName;
  std::vector<std::string> itsDirections;

  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
