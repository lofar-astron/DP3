// Counter.cc: DP3 step class to count flags
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Counter.h"

#include <iostream>
#include <fstream>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Counter::Counter(const common::ParameterSet& parset, const string& prefix)
    : name_(prefix),
      count_(0),
      save_to_json_(parset.getBool(prefix + "savetojson", false)),
      json_filename_(parset.getString(prefix + "jsonfilename",
                                      "FlagPercentagePerStation.JSON")),
      flag_counter_(parset, prefix) {}

Counter::~Counter() {}

void Counter::show(std::ostream& os) const {
  os << "Counter " << name_ << '\n';
}

void Counter::showCounts(std::ostream& os) const {
  os << "\nCumulative flag counts in Counter " << name_;
  os << "\n=================================\n";
  flag_counter_.showBaseline(os, count_);
  flag_counter_.showChannel(os, count_);
  if (save_to_json_) {
    os << "\nSaving counts to JSON file " << json_filename_ << "\n";

    std::ostringstream percentage_per_station;
    flag_counter_.showStation(percentage_per_station, count_);
    std::ofstream jsonfile(json_filename_);
    jsonfile << percentage_per_station.str();
    jsonfile.close();
  }
}

void Counter::updateInfo(const base::DPInfo& info_in) {
  Step::updateInfo(info_in);
  // Initialize the flag counters.
  flag_counter_.init(info_in);
}

bool Counter::process(const base::DPBuffer& buf) {
  const casacore::IPosition& shape = buf.GetCasacoreFlags().shape();
  unsigned int nrcorr = shape[0];
  unsigned int nrchan = shape[1];
  unsigned int nrbl = shape[2];
  const bool* flag_pointer = buf.GetFlags().data();
  for (unsigned int i = 0; i < nrbl; ++i) {
    for (unsigned int j = 0; j < nrchan; ++j) {
      if (*flag_pointer) {
        flag_counter_.incrBaseline(i);
        flag_counter_.incrChannel(j);
      }
      flag_pointer += nrcorr;  // only count 1st corr
    }
  }
  // Let the next step do its processing.
  getNextStep()->process(buf);
  ++count_;
  return true;
}

void Counter::finish() {
  // Let the next step finish its processing.
  getNextStep()->finish();
}

}  // namespace steps
}  // namespace dp3
