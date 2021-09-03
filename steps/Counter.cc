// Counter.cc: DPPP step class to count flags
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Counter.h"

#include "../base/DPInfo.h"

#include "../common/ParameterSet.h"

#include <iostream>
#include <fstream>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Counter::Counter(InputStep* input, const common::ParameterSet& parset,
                 const string& prefix)
    : itsName(prefix),
      itsCount(0),
      itsSaveToJson(parset.getBool(prefix + "savetojson", false)),
      itsJsonFilename(parset.getString(prefix + "jsonfilename",
                                       "FlagPercentagePerStation.JSON")),
      itsFlagCounter(input->msName(), parset, prefix) {
  itsFlagData = parset.getBool(prefix + "flagdata", false);
}

Counter::~Counter() {}

void Counter::show(std::ostream& os) const {
  os << "Counter " << itsName << '\n';
}

void Counter::showCounts(std::ostream& os) const {
  os << "\nCumulative flag counts in Counter " << itsName;
  os << "\n=================================\n";
  itsFlagCounter.showBaseline(os, itsCount);
  itsFlagCounter.showChannel(os, itsCount);
  if (itsSaveToJson) {
    os << "\nSaving counts to JSON file " << itsJsonFilename << "\n";

    std::ostringstream percentagePerStation;
    itsFlagCounter.showStation(percentagePerStation, itsCount);
    std::ofstream jsonfile(itsJsonFilename);
    jsonfile << percentagePerStation.str();
    jsonfile.close();
  }
}

void Counter::updateInfo(const base::DPInfo& infoIn) {
  info() = infoIn;
  // Visibility data must be read if needed, so NaNs are flagged.
  if (itsFlagData) {
    info().setNeedVisData();
  }
  // Initialize the flag counters.
  itsFlagCounter.init(infoIn);
}

bool Counter::process(const base::DPBuffer& buf) {
  const casacore::IPosition& shape = buf.getFlags().shape();
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

}  // namespace steps
}  // namespace dp3
