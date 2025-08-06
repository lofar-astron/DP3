// Step.cc: Abstract base class for a DP3 step
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/steps/Step.h>

#include <cassert>

#include <aocommon/system.h>
#include <aocommon/threadpool.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;

namespace dp3 {
namespace steps {

Step::Step() {
  if (!threading_is_initialized_) {
    aocommon::ThreadPool::GetInstance().SetNThreads(
        aocommon::system::ProcessorCount());
    threading_is_initialized_ = true;
  }
}

Step::~Step() = default;

void Step::setInfo(const DPInfo& info) {
  // Update the info of this step using the given info.
  updateInfo(info);
  // If there is a next step, set its info using the info of this step.
  if (getNextStep()) {
    getNextStep()->setInfo(getInfoOut());
  }
}

void Step::updateInfo(const DPInfo& info) {
  input_info_ = info;
  output_info_ = info;
}

void Step::addToMS(const std::string& msName) {
  if (previous_step_) previous_step_->addToMS(msName);
}

void Step::showCounts(std::ostream&) const {}

void Step::showTimings(std::ostream&, double) const {}

}  // namespace steps
}  // namespace dp3
