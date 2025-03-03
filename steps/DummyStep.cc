// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DummyStep.h"

#include <iostream>

#include "../base/FlagCounter.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

DummyStep::DummyStep([[maybe_unused]] const common::ParameterSet& parset,
                     const std::string& prefix)
    : name_(prefix) {}

void DummyStep::updateInfo(const DPInfo& info_in) { Step::updateInfo(info_in); }

void DummyStep::show(std::ostream& os) const {
  os << "DummyStep " << name_ << '\n';
}

void DummyStep::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " DummyStep " << name_ << '\n';
}

bool DummyStep::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  // Do the regular processing here. 'buffer' will contain at least
  // the required fields of the step, as returned from getRequiredFields.

  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

bool DummyStep::process(std::unique_ptr<base::BdaBuffer> buffer) {
  // Do the BDA processing here.
  return getNextStep()->process(std::move(buffer));
}

void DummyStep::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
