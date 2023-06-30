// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_NULLSTEP_H_
#define DP3_STEPS_NULLSTEP_H_

#include "OutputStep.h"

namespace dp3 {
namespace steps {

/// @brief This class defines a null step in the DP3 pipeline.
/// It can be used as the last step in the pipeline, so other steps
/// do not need to test if there is a next step.
class NullStep : public OutputStep {
 public:
  ~NullStep() override {}

  /// A null step requires nothing.
  common::Fields getRequiredFields() const override { return {}; }

  /// A null step provides nothing.
  common::Fields getProvidedFields() const override { return {}; }

  /// Process regular data. It does nothing.
  bool process(std::unique_ptr<base::DPBuffer>) override { return true; }

  /// Process bda data. It does nothing.
  bool process(std::unique_ptr<base::BDABuffer>) override { return true; }

  /// Finish the processing of this step and subsequent steps.
  /// It does nothing.
  void finish() override {}

  /// Show the step parameters.
  /// It does nothing.
  void show(std::ostream&) const override {}

  /// Accept BDA and regular data.
  bool accepts(MsType t) const override {
    return t == MsType::kRegular || t == MsType::kBda;
  }
};

}  // namespace steps
}  // namespace dp3

#endif
