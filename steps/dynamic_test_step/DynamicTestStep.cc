// TestDyn.cc: Test of a dynamically loaded DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

// @file
// @brief Test of a dynamically loaded DPPP step

#include "DynamicTestStep.h"

#include "../../base/DPBuffer.h"
#include "../../base/DP3.h"

namespace dp3 {
namespace steps {
namespace dynamic_test_step {

DynamicTestStep::DynamicTestStep(InputStep* input,
                                 const common::ParameterSet& pset,
                                 const std::string& prefix)
    : Averager(*input, pset, prefix) {}

DynamicTestStep::~DynamicTestStep() {}

Step::ShPtr DynamicTestStep::makeStep(InputStep* input,
                                      const common::ParameterSet& pset,
                                      const std::string& prefix) {
  return std::make_shared<DynamicTestStep>(input, pset, prefix);
}

}  // namespace dynamic_test_step
}  // namespace steps
}  // namespace dp3

// Define the function to make the TestDynStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
void register_testdyndppp() {
  dp3::base::DP3::registerStepCtor(
      "TestDynDPPP", dp3::steps::dynamic_test_step::DynamicTestStep::makeStep);
}
