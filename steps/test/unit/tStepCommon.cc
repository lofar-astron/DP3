// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper routines for testing steps.

#include "tStepCommon.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"

namespace dp3 {
namespace steps {
namespace test {

/// Connect a series of steps and execute them.
void Execute(const std::vector<std::shared_ptr<Step>>& steps) {
  // Connect the steps.
  for (std::size_t i = 1; i < steps.size(); ++i) {
    steps[i - 1]->setNextStep(steps[i]);
  }

  // Set DPInfo for all steps.
  steps.front()->setInfo(base::DPInfo());

  // Finally, execute the steps.
  base::DPBuffer buf;
  while (steps.front()->process(buf))
    ;
  steps.front()->finish();
}

dp3::common::ParameterSet CreateParameterSet(
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  dp3::common::ParameterSet result;
  for (const auto& parameter : parameters)
    result.add(parameter.first, parameter.second);
  return result;
}

}  // namespace test
}  // namespace steps
}  // namespace dp3
