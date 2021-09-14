// GainCal.cc: DPPP step class to DummyStep visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "DummyStep.h"

#include <iostream>

#include "../base/BDABuffer.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

DummyStep::DummyStep(InputStep* input,
                     [[maybe_unused]] const common::ParameterSet& parset,
                     [[maybe_unused]] const string& prefix)
    : itsInput(input) {}

DummyStep::~DummyStep() {}

void DummyStep::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
}

void DummyStep::show(std::ostream& os) const {
  os << "DummyStep " << itsName << '\n';
}

void DummyStep::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " DummyStep " << itsName << '\n';
}

bool DummyStep::process(const DPBuffer& bufin) {
  itsTimer.start();
  itsBuffer.copy(bufin);
  itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
  itsInput->fetchWeights(bufin, itsBuffer, itsTimer);

  itsTimer.stop();
  getNextStep()->process(itsBuffer);
  return false;
}

bool DummyStep::process(std::unique_ptr<base::BDABuffer> buffer) {
  return getNextStep()->process(std::move(buffer));
}

void DummyStep::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
