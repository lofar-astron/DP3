// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ResultStep.h"

#include "NullStep.h"

namespace dp3 {
namespace steps {

ResultStep::ResultStep() : buffer_() {
  setNextStep(std::make_shared<NullStep>());
}

}  // namespace steps
}  // namespace dp3