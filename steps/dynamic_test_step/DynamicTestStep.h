// TestDynStep.h: Test of a dynamically loaded DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Test of a dynamically loaded DPPP step
/// @author Ger van Diepen

#ifndef TESTDYNDPPP_TESTDYNSTEP_H
#define TESTDYNDPPP_TESTDYNSTEP_H

#include "../../steps/Step.h"
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
  DynamicTestStep(InputStep*, const common::ParameterSet&, const std::string&);
  virtual ~DynamicTestStep();
  static Step::ShPtr makeStep(InputStep*, const common::ParameterSet&,
                              const std::string&);
};

}  // namespace dynamic_test_step
}  // namespace steps
}  // namespace dp3

// Define the function (without name mangling) to register the 'constructor'.
extern "C" {
void register_testdyndppp();
}

#endif
