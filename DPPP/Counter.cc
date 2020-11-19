// Counter.cc: DPPP step class to count flags
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Counter.h"
#include "DPInfo.h"

#include "../Common/ParameterSet.h"

#include <iostream>

using namespace casacore;

namespace DP3 {
namespace DPPP {

Counter::Counter(DPInput* input, const ParameterSet& parset,
                 const string& prefix)
    : itsName(prefix),
      itsCount(0),
      itsFlagCounter(input->msName(), parset, prefix) {
  itsFlagData = parset.getBool(prefix + "flagdata", false);
}

Counter::~Counter() {}

void Counter::show(std::ostream& os) const {
  os << "Counter " << itsName << std::endl;
}

void Counter::showCounts(std::ostream& os) const {
  os << endl << "Cumulative flag counts in Counter " << itsName;
  os << endl << "=================================" << endl;
  itsFlagCounter.showBaseline(os, itsCount);
  itsFlagCounter.showChannel(os, itsCount);
}

void Counter::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  // Visibility data must be read if needed, so NaNs are flagged.
  if (itsFlagData) {
    info().setNeedVisData();
  }
  // Initialize the flag counters.
  itsFlagCounter.init(infoIn);
}

bool Counter::process(const DPBuffer& buf) {
  const IPosition& shape = buf.getFlags().shape();
  unsigned int nrcorr = shape[0];
  unsigned int nrchan = shape[1];
  unsigned int nrbl = shape[2];
  const bool* flagPtr = buf.getFlags().data();
  for (unsigned int i = 0; i < nrbl; ++i) {
    for (unsigned int j = 0; j < nrchan; ++j) {
      if (*flagPtr) {
        itsFlagCounter.incrBaseline(i);
        itsFlagCounter.incrChannel(j);
      }
      flagPtr += nrcorr;  // only count 1st corr
    }
  }
  // Let the next step do its processing.
  getNextStep()->process(buf);
  itsCount++;
  return true;
}

void Counter::finish() {
  // Let the next step finish its processing.
  getNextStep()->finish();
}

}  // namespace DPPP
}  // namespace DP3
