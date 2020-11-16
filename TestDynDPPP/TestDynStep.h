// TestDynStep.h: Test of a dynamically loaded DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Test of a dynamically loaded DPPP step
/// @author Ger van Diepen

#ifndef TESTDYNDPPP_TESTDYNSTEP_H
#define TESTDYNDPPP_TESTDYNSTEP_H

#include "../DPPP/DPStep.h"
#include "../DPPP/Averager.h"
#include "../DPPP/DPInput.h"
#include "../Common/ParameterSet.h"

namespace DP3 {
namespace DPPP {
/// @brief Test of a dynamically loaded DPPP step

/// This class is a test (and an example) of a DPStep loaded
/// dynamically from a shared library.
/// To make test life easy it uses the Averager class underneath.

class TestDynStep : public Averager {
 public:
  TestDynStep(DPInput*, const ParameterSet&, const std::string&);
  virtual ~TestDynStep();
  static DPStep::ShPtr makeStep(DPInput*, const ParameterSet&,
                                const std::string&);
};

}  // namespace DPPP
}  // namespace DP3

// Define the function (without name mangling) to register the 'constructor'.
extern "C" {
void register_testdyndppp();
}

#endif
