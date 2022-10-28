// TestDyn.cc: Test of a dynamically loaded DP3 step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

// @file
// @brief Test of a dynamically loaded DP3 step

#include "DynamicTestStep.h"

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DP3.h>

namespace dp3 {
namespace steps {
namespace dynamic_test_step {

DynamicTestStep::DynamicTestStep(const common::ParameterSet& parset,
                                 const std::string& prefix)
    : Averager(parset, prefix) {}

std::shared_ptr<Step> DynamicTestStep::MakeStep(
    const common::ParameterSet& parset, const std::string& prefix) {
  return std::make_shared<DynamicTestStep>(parset, prefix);
}

}  // namespace dynamic_test_step
}  // namespace steps
}  // namespace dp3

// Define the function to make the TestDynStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
void register_testdyndp3() {
  dp3::base::DP3::RegisterStepCreator(
      "TestDynDP3", dp3::steps::dynamic_test_step::DynamicTestStep::MakeStep);
}
