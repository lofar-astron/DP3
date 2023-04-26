// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionInterval.h"

namespace dp3 {
namespace base {

SolutionInterval::SolutionInterval(const std::size_t buffer_size)
    : buffers_(), original_flags_(buffer_size), original_weights_(buffer_size) {
  buffers_.reserve(buffer_size);
}

void SolutionInterval::PushBack(std::unique_ptr<DPBuffer> buffer) {
  if (buffers_.size() >= original_flags_.size()) {
    throw std::runtime_error("SolutionInterval exceeds buffer size");
  }

  original_flags_[buffers_.size()].assign(buffer->GetCasacoreFlags());
  original_weights_[buffers_.size()].assign(buffer->GetCasacoreWeights());
  // Ensure that the copies are independent of the data in 'buffer'.
  original_flags_[buffers_.size()].unique();
  original_weights_[buffers_.size()].unique();

  buffers_.push_back(std::move(buffer));
}

void SolutionInterval::RestoreFlagsAndWeights() {
  constexpr common::Fields kFlagsField(common::Fields::Single::kFlags);
  constexpr common::Fields kWeightsField(common::Fields::Single::kWeights);
  const common::Fields kFields(kFlagsField | kWeightsField);
  for (std::size_t index = 0; index < buffers_.size(); ++index) {
    buffers_[index]->GetCasacoreFlags() = original_flags_[index];
    buffers_[index]->GetCasacoreWeights() = original_weights_[index];
    buffers_[index]->MakeIndependent(kFields);
  }
}

}  // namespace base
}  // namespace dp3
