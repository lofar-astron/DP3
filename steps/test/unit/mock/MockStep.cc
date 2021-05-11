// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockStep.h"

#include "../../../../base/BDABuffer.h"

#include <boost/test/unit_test.hpp>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;

namespace dp3 {
namespace steps {

MockStep::MockStep(std::function<void(const DPBuffer&)>* check_buffer)
    : check_buffer_(check_buffer), bda_buffers_(), finish_count_(0) {}

MockStep::~MockStep() {}

bool MockStep::process(const DPBuffer& buffer) {
  // Ask TJ: if the BOOST_CHECK is present, any call to BDAExpander->GetNextStep
  // process fails (during unit tests) BOOST_CHECK(check_buffer_);
  // (*check_buffer_)(buffer);
  regular_buffers_.push_back(std::move(buffer));
  return true;
}

bool MockStep::process(std::unique_ptr<BDABuffer> buffer) {
  bda_buffers_.push_back(std::move(buffer));
  return true;
}

void MockStep::show(std::ostream&) const { BOOST_CHECK(false); }

void MockStep::ClearBdaBuffers() { bda_buffers_.clear(); }

std::size_t MockStep::TotalRowCount() const {
  std::size_t count = 0;
  for (const std::unique_ptr<BDABuffer>& buffer : bda_buffers_) {
    count += buffer->GetRows().size();
  }
  return count;
}

}  // namespace steps
}  // namespace dp3
