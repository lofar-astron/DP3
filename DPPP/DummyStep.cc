// GainCal.cc: DPPP step class to DummyStep visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "DummyStep.h"

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace DP3 {
namespace DPPP {

DummyStep::DummyStep(DPInput* input, const ParameterSet& parset,
                     const string& prefix)
    : itsInput(input) {}

DummyStep::~DummyStep() {}

void DummyStep::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
}

void DummyStep::show(std::ostream& os) const {
  os << "DummyStep " << itsName << endl;
}

void DummyStep::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " DummyStep " << itsName << endl;
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

void DummyStep::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace DPPP
}  // namespace DP3
