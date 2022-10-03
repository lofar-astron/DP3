// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MultiResultStep.h"

#include "NullStep.h"

namespace dp3 {
namespace steps {

MultiResultStep::MultiResultStep(unsigned int size) : buffers_(size), size_(0) {
  setNextStep(std::make_shared<NullStep>());
}

bool MultiResultStep::process(const base::DPBuffer& buffer) {
  buffers_[size_].copy(buffer);
  ++size_;
  getNextStep()->process(buffer);
  return true;
}

}  // namespace steps
}  // namespace dp3