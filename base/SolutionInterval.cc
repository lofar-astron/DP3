// SolutionInterval.cc
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionInterval.h"

using dp3::steps::InputStep;

namespace dp3 {
namespace base {

SolutionInterval::SolutionInterval(InputStep& input,
                                   const std::size_t n_solution,
                                   const std::size_t buffer_size,
                                   common::NSTimer timer)
    : buffer_size_(buffer_size),
      n_solution_(n_solution),
      timer_(timer),
      input_(input),
      buffer_index_(0),
      buffers_(buffer_size),
      original_flags_(buffer_size),
      original_weights_(buffer_size) {}

SolutionInterval::~SolutionInterval() {}

void SolutionInterval::PushBack(const DPBuffer& buffer) {
  if (buffer_index_ >= buffers_.size()) {
    throw std::runtime_error("SolutionInterval exceeds buffer size");
  }

  input_.fetchUVW(buffer, buffers_[buffer_index_], timer_);
  input_.fetchWeights(buffer, buffers_[buffer_index_], timer_);
  input_.fetchFullResFlags(buffer, buffers_[buffer_index_], timer_);

  buffers_[buffer_index_].copy(buffer);

  original_flags_[buffer_index_].assign(buffer.getFlags());
  original_weights_[buffer_index_].assign(buffer.getWeights());

  ++buffer_index_;
}

void SolutionInterval::RestoreFlagsAndWeights() {
  for (std::size_t index = 0; index < buffer_index_; ++index) {
    buffers_[index].getFlags().assign(original_flags_[index]);
    buffers_[index].getWeights().assign(original_weights_[index]);
  }
}

}  // namespace base
}  // namespace dp3
