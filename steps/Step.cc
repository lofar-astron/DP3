// Step.cc: Abstract base class for a DP3 step
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/steps/Step.h>

#include <assert.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;

namespace dp3 {
namespace steps {

Step::~Step() {}

const DPInfo& Step::setInfo(const DPInfo& info) {
  // Update the info of this step using the given info.
  updateInfo(info);
  // If there is a next step, set its info using the info of this step.
  if (getNextStep()) {
    return getNextStep()->setInfo(getInfo());
  }
  return getInfo();
}

void Step::updateInfo(const DPInfo& infoIn) { info() = infoIn; }

void Step::addToMS(const string& msName) {
  if (itsPrevStep) itsPrevStep->addToMS(msName);
}

void Step::showCounts(std::ostream&) const {}

void Step::showTimings(std::ostream&, double) const {}

}  // namespace steps
}  // namespace dp3
