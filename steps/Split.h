// Split.h: DP3 step class to Split visibilities from a source model
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to Split visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_SPLIT_H_
#define DP3_STEPS_SPLIT_H_

#include <utility>

#include <dp3/base/DPBuffer.h>

#include "../common/ParameterSet.h"

#include "InputStep.h"
#include "OutputStep.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step class to Split visibilities from a source model
class Split : public OutputStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Split(InputStep*, const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override;

  void SetFieldsToWrite(const common::Fields& fields) override;

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Override setNextStep, since Split should be the last step.
  void setNextStep(std::shared_ptr<Step> step) override;

 private:
  std::string name_;

  /// The names of the parameters that differ along the substeps.
  std::vector<std::string> replace_parameters_;

  /// The first step in each chain of sub steps.
  std::vector<std::shared_ptr<Step>> sub_steps_;
};

}  // namespace steps
}  // namespace dp3

#endif
