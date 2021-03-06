// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper routines for testing steps.

#include "../../Step.h"
#include <vector>

namespace dp3 {
namespace steps {
namespace test {

/// Connect a series of steps and execute them.
void Execute(const std::vector<std::shared_ptr<Step>>& steps);

}  // namespace test
}  // namespace steps
}  // namespace dp3