// Step.cc: Abstract base class for a DP3 step
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "steps/Step.h"

#include <cassert>

#include <aocommon/system.h>

#include <schaapcommon/threading/threadpool.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;

namespace dp3 {
namespace steps {

Step::Step() {
  if (!threading_is_initialized_) {
    schaapcommon::ThreadPool::GetInstance().SetNThreads(
        aocommon::system::ProcessorCount());
    threading_is_initialized_ = true;
  }
}

Step::~Step() {
  // If the step has a previous step, that previous step's next step pointer
  // should point to this step, which means there is still a valid pointer
  // to this step and this destructor can't be called.
  assert(!previous_step_);

  // If the step has a next step, update that next step.
  if (next_step_) {
    assert(next_step_->previous_step_ == this);
    next_step_->previous_step_ = nullptr;
  }
}

void Step::setNextStep(std::shared_ptr<Step> next_step) {
  // Update any existing next step, so it can be safely destroyed, without
  // triggering any assertion in the Step destructor.
  if (next_step_) {
    assert(next_step_->previous_step_ == this);
    next_step_->previous_step_ = nullptr;
  }

  next_step_ = std::move(next_step);
  if (next_step_) {
    next_step_->previous_step_ = this;
  }
}

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
