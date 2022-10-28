// TestDynStep.h: Test of a dynamically loaded DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Test of a dynamically loaded DPPP step
/// @author Ger van Diepen

#ifndef DP3_STEPS_DYNAMICTESTSTEP_H_
#define DP3_STEPS_DYNAMICTESTSTEP_H_

#include <dp3/steps/Step.h>
#include "../../steps/Averager.h"
#include "../../steps/InputStep.h"
#include "../../common/ParameterSet.h"

namespace dp3 {
namespace steps {
namespace dynamic_test_step {
/// @brief Test of a dynamically loaded DPPP step

/// This class is a test (and an example) of a Step loaded
/// dynamically from a shared library.
/// To make test life easy it uses the Averager class underneath.

class DynamicTestStep : public Averager {
 public:
  DynamicTestStep(const common::ParameterSet&, const std::string&);
  static std::shared_ptr<Step> MakeStep(const common::ParameterSet&,
                                        const std::string&);
};

}  // namespace dynamic_test_step
}  // namespace steps
}  // namespace dp3

// Define the function (without name mangling) to register the 'constructor'.
extern "C" {
void register_testdyndp3();
}

#endif
