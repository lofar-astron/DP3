// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../../base/Direction.h"

#include <string>
#include <utility>

namespace dp3 {
namespace steps {
namespace test {

/// MS name for the predict tests.
const std::string kPredictSourceDB = "tNDPPP-generic.MS/sky";

/// Expected first direction when using tNDPPP-generic.MS.
/// Multiple predict tests use this value.
const base::Direction kExpectedFirstDirection{0.0097549942552467052,
                                              0.55260038817991308};
}  // namespace test
}  // namespace steps
}  // namespace dp3
