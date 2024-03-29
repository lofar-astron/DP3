// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper routines for testing steps.

#include "tStepCommon.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include <aocommon/threadpool.h>

namespace dp3 {
namespace steps {
namespace test {

/// Connect a series of steps and execute them.
void Execute(const std::vector<std::shared_ptr<Step>>& steps) {
  // Connect the steps.
  for (std::size_t i = 1; i < steps.size(); ++i) {
    steps[i - 1]->setNextStep(steps[i]);
  }

  // Set DPInfo for all steps. Use a single thread in tests.
  base::DPInfo info;
  aocommon::ThreadPool::GetInstance().SetNThreads(1);
  steps.front()->setInfo(info);

  // Finally, execute the steps.
  while (steps.front()->process(std::make_unique<base::DPBuffer>()))
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
