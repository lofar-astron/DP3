// SolutionInterval.cc
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionInterval.h"

using dp3::steps::InputStep;

namespace dp3 {
namespace base {

SolutionInterval::SolutionInterval(const std::size_t n_solution,
                                   const std::size_t buffer_size,
                                   common::NSTimer timer)
    : buffer_size_(buffer_size),
      n_solution_(n_solution),
      timer_(timer),
      buffer_index_(0),
      buffers_(buffer_size),
      original_flags_(buffer_size),
      original_weights_(buffer_size) {}

SolutionInterval::~SolutionInterval() {}

void SolutionInterval::PushBack(const DPBuffer& buffer) {
  if (buffer_index_ >= buffers_.size()) {
    throw std::runtime_error("SolutionInterval exceeds buffer size");
  }

  buffers_[buffer_index_].copy(buffer);

  original_flags_[buffer_index_].assign(buffer.GetCasacoreFlags());
  original_weights_[buffer_index_].assign(buffer.GetCasacoreWeights());
  // Ensure that the copies are independent of the data in 'buffer'..
  original_flags_[buffer_index_].unique();
  original_weights_[buffer_index_].unique();

  ++buffer_index_;
}

void SolutionInterval::RestoreFlagsAndWeights() {
  constexpr common::Fields kFlagsField(common::Fields::Single::kFlags);
  constexpr common::Fields kWeightsField(common::Fields::Single::kWeights);
  const common::Fields kFields(kFlagsField | kWeightsField);
  for (std::size_t index = 0; index < buffer_index_; ++index) {
    buffers_[index].GetCasacoreFlags() = original_flags_[index];
    buffers_[index].GetCasacoreWeights() = original_weights_[index];
    buffers_[index].MakeIndependent(kFields);
  }
}

}  // namespace base
}  // namespace dp3
