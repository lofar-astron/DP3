// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_TEST_UNIT_H5PARMFIXTURE_H_
#define DP3_STEPS_TEST_UNIT_H5PARMFIXTURE_H_

#include <filesystem>

#include <schaapcommon/h5parm/h5parm.h>

#include "tPredict.h"

namespace dp3 {
namespace steps {
namespace test {

/// Creates a basic H5Parm file for testing H5ParmPredict.
struct H5ParmFixture {
  H5ParmFixture() {
    schaapcommon::h5parm::H5Parm h5parm(kParmDb, true);
    // "gain" must be accepted by JonesParameters::StringToCorrectType.
    // "dir" is the name of a required h5parm axis.
    auto soltab = h5parm.CreateSolTab(kSoltabName, "gain", {{"dir", 1}});
    soltab.SetSources({'[' + kPredictDirection + ']'});
    // Set a dummy value and weight.
    soltab.SetValues({42.0}, {1.0}, "CREATE with H5ParmFixture");
  }

  ~H5ParmFixture() { std::filesystem::remove(kParmDb); }

  const std::string kParmDb = "test_h5parmfixture.h5";
  const std::string kSoltabName = "foobar42";
};

}  // namespace test
}  // namespace steps
}  // namespace dp3

#endif