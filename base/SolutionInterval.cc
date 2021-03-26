// SolutionInterval.cc
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionInterval.h"

using dp3::steps::InputStep;

namespace dp3 {
namespace base {

SolutionInterval::SolutionInterval(InputStep* input,
                                   const std::size_t n_solution,
                                   const std::size_t buffer_size,
                                   const std::size_t n_dir,
                                   common::NSTimer timer)
    : buffer_size_(buffer_size),
      n_solution_(n_solution),
      timer_(timer),
      input_(input),
      buffer_index_(0),
      buffers_(buffer_size),
      data_ptrs_(buffer_size),
      weight_ptrs_(buffer_size),
      original_flags_(buffer_size),
      original_weights_(buffer_size),
      model_data_ptrs_(buffer_size),
      model_data_(buffer_size) {
  for (std::size_t t = 0; t < buffer_size_; ++t) {
    model_data_ptrs_[t].resize(n_dir);
  }
}

SolutionInterval::~SolutionInterval() {}

void SolutionInterval::CopyBuffer(const DPBuffer& buffer) {
  if (buffer_index_ >= buffers_.capacity()) {
    throw std::runtime_error("SolutionInterval exceeds buffer size");
  }

  input_->fetchUVW(buffer, buffers_[buffer_index_], timer_);
  input_->fetchWeights(buffer, buffers_[buffer_index_], timer_);
  input_->fetchFullResFlags(buffer, buffers_[buffer_index_], timer_);

  buffers_[buffer_index_].copy(buffer);

  original_flags_[buffer_index_].assign(buffer.getFlags());
  original_weights_[buffer_index_].assign(buffer.getWeights());
  data_ptrs_[buffer_index_] = buffers_[buffer_index_].getData().data();
  weight_ptrs_[buffer_index_] = buffers_[buffer_index_].getWeights().data();

  ++buffer_index_;
}

void SolutionInterval::RestoreFlagsAndWeights() {
  for (std::size_t index = 0; index < buffer_index_; ++index) {
    buffers_[index].getFlags().assign(original_flags_[index]);
    buffers_[index].getWeights().assign(original_weights_[index]);
  }
}

void SolutionInterval::Fit() {
  data_ptrs_.resize(buffer_index_);
  weight_ptrs_.resize(buffer_index_);
  model_data_.resize(buffer_index_);
  model_data_ptrs_.resize(buffer_index_);
}
}  // namespace base
}  // namespace dp3
