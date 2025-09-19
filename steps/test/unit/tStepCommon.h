// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper routines for testing steps.

#ifndef DP3_STEPS_TEST_UNIT_TSTEPCOMMON_H_
#define DP3_STEPS_TEST_UNIT_TSTEPCOMMON_H_

#include <memory>
#include <vector>

#include <dp3/steps/Step.h>
#include "../../../common/ParameterSet.h"

namespace dp3 {
namespace steps {
namespace test {

/// Connect a series of steps and execute them.
void Execute(const std::vector<std::shared_ptr<Step>>& steps);

/** Returns a @ref dp3::common::ParameterSet based on @a parameters. */
dp3::common::ParameterSet CreateParameterSet(
    const std::vector<std::pair<std::string, std::string>>& parameters);

/**
 * Helper function to aid testing @c Step::show.
 * @return The result of Step::show, using the C locale.
 */
inline std::string Show(const Step& step) {
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
