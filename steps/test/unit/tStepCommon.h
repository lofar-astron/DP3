// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper routines for testing steps.

#ifndef DPPP_STEPS_TEST_UNIT_TSTEPCOMMON_H
#define DPPP_STEPS_TEST_UNIT_TSTEPCOMMON_H

#include "../../Step.h"
#include "../../InputStep.h"
#include "../../../common/ParameterSet.h"
#include <vector>

namespace dp3 {
namespace steps {
namespace test {

/// Connect a series of steps and execute them.
void Execute(const std::vector<std::shared_ptr<Step>>& steps);

/**
 * Helper function to aid testing @c Step::show.
 *
 * Creates a @a Step using @a parameters for its parset and a fixed prefix
 * @c "prefix.". It returns the result of @a Step::show.
 */
template <class Step>
inline std::string Show(
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  class : public dp3::steps::InputStep {
    void finish() override {}
    void show(std::ostream&) const override {}
  } input;

  dp3::common::ParameterSet parset;
  for (const auto& parameter : parameters)
    parset.add(parameter.first, parameter.second);

  const Step step{&input, parset, "prefix."};
  std::stringstream output;
  // Ensure the test doesn't depend on the system's locale settings.
  output.imbue(std::locale::classic());
  step.show(output);
  return output.str();
}

}  // namespace test
}  // namespace steps
}  // namespace dp3

#endif
