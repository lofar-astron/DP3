// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_TEST_UNIT_TPREDICT_H_
#define DP3_STEPS_TEST_UNIT_TPREDICT_H_

#include <string>

#include "base/Direction.h"
#include "test_config.h"

namespace dp3 {
namespace steps {
namespace test {

/// MS name for the predict tests.
const std::string kPredictSkymodel =
    DP3_RESOURCE_DIR "/tNDPPP-generic-skymodel.txt";
const std::string kPredictDirection = "0002.2+3139";

/// Expected first direction when using tNDPPP-generic.MS.
/// Multiple predict tests use this value.
constexpr base::Direction kExpectedFirstDirection{0.0097549942552467052,
                                                  0.55260038817991308};
}  // namespace test
}  // namespace steps
}  // namespace dp3

#endif
