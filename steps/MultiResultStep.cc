// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MultiResultStep.h"

#include <cassert>

#include "NullStep.h"

namespace dp3 {
namespace steps {

MultiResultStep::MultiResultStep(unsigned int size) : buffers_(size), size_(0) {
  setNextStep(std::make_shared<NullStep>());
}

bool MultiResultStep::process(std::unique_ptr<base::DPBuffer> buffer) {
  assert(size_ < buffers_.size());

  // If the next step is not a NullStep, one copy of the buffer is saved into
  // buffers_[size_] and the other is moved to the next step's process()
  // function. Remove the MakeIndependent() function once the casacore cubes are
  // replaced by xtensor.
  if (dynamic_cast<NullStep*>(getNextStep().get()) == nullptr) {
    buffers_[size_] = std::make_unique<base::DPBuffer>(*buffer);
    buffers_[size_]->MakeIndependent(kDataField | kWeightsField | kFlagsField |
                                     kUvwField);
    ++size_;
    getNextStep()->process(std::move(buffer));
  } else {
    buffers_[size_] = std::move(buffer);
    buffers_[size_]->MakeIndependent(kDataField | kWeightsField | kFlagsField |
                                     kUvwField);
    ++size_;
  }
  return true;
}

}  // namespace steps
}  // namespace dp3