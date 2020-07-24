// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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

#include "MockStep.h"

#include "../../../BDABuffer.h"

#include <boost/test/unit_test.hpp>

namespace DP3 {
namespace DPPP {

MockStep::MockStep(std::function<void(const DPBuffer&)>* check_buffer)
    : check_buffer_(check_buffer), bda_buffers_(), finish_count_(0) {}

MockStep::~MockStep() {}

bool MockStep::process(const DPBuffer& buffer) {
  BOOST_CHECK(check_buffer_);
  (*check_buffer_)(buffer);
  return true;
}

bool MockStep::process(std::unique_ptr<BDABuffer> buffer) {
  bda_buffers_.push_back(std::move(buffer));
  return true;
}

void MockStep::show(std::ostream&) const { BOOST_CHECK(false); }

void MockStep::ClearBdaBuffers() { bda_buffers_.clear(); }

void MockStep::CheckFinishCount(std::size_t expected_count) const {
  BOOST_CHECK_EQUAL(finish_count_, expected_count);
}

}  // namespace DPPP
}  // namespace DP3
