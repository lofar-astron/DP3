// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockStep.h"

#include "base/BdaBuffer.h"

#include <boost/test/unit_test.hpp>

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;

namespace dp3 {
namespace steps {

std::size_t MockStep::TotalRowCount() const {
  std::size_t count = 0;
  for (const std::unique_ptr<BdaBuffer>& buffer : bda_buffers_) {
    count += buffer->GetRows().size();
  }
  return count;
}

}  // namespace steps
}  // namespace dp3
