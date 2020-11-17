// SolutionInterval.cc
//
// Copyright (C) 2020
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "SolutionInterval.h"

namespace DP3 {

namespace DPPP {

SolutionInterval::SolutionInterval(DPInput* input, const std::size_t n_solution,
                                   const std::size_t buffer_size,
                                   const std::size_t n_dir, NSTimer timer)
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
}  // namespace DPPP
}  // namespace DP3
