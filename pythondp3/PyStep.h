// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_PYDPSTEP_H_
#define DP3_STEPS_PYDPSTEP_H_

#include "steps/Step.h"
#include "../common/ParameterSet.h"

namespace dp3 {
namespace pythondp3 {

class PyStep final : public steps::Step {
 public:
  static std::shared_ptr<PyStep> create_instance(
      const common::ParameterSet& parset, const std::string& prefix);
  using steps::Step::Step;

  PyStep();

  void show(std::ostream& os) const override;

  common::Fields getRequiredFields() const override;
  common::Fields getProvidedFields() const override;

  void updateInfo(const base::DPInfo&) override;

  bool process(std::unique_ptr<base::DPBuffer>) override;

  void finish() override;
};

}  // namespace pythondp3
}  // namespace dp3

#endif  // DP3_STEPS_PYDPSTEP_H_
