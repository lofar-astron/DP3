// TestDyn.cc: Test of a dynamically loaded DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

// @file
// @brief Test of a dynamically loaded DPPP step

#include "TestDynStep.h"

#include "../DPPP/DPBuffer.h"
#include "../DPPP/DPRun.h"

namespace DP3 {
namespace DPPP {

TestDynStep::TestDynStep(DPInput* input, const ParameterSet& pset,
                         const std::string& prefix)
    : Averager(input, pset, prefix) {}

TestDynStep::~TestDynStep() {}

DPStep::ShPtr TestDynStep::makeStep(DPInput* input, const ParameterSet& pset,
                                    const std::string& prefix) {
  return std::make_shared<TestDynStep>(input, pset, prefix);
}

}  // namespace DPPP
}  // namespace DP3

// Define the function to make the TestDynStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
void register_testdyndppp() {
  DP3::DPPP::DPRun::registerStepCtor("TestDynDPPP",
                                     DP3::DPPP::TestDynStep::makeStep);
}
