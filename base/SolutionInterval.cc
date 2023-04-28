// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionInterval.h"

namespace dp3 {
namespace base {

SolutionInterval::SolutionInterval(const std::size_t buffer_size) : buffers_() {
  buffers_.reserve(buffer_size);
}

void SolutionInterval::PushBack(std::unique_ptr<DPBuffer> buffer) {
  buffers_.push_back(std::move(buffer));
}

}  // namespace base
}  // namespace dp3
