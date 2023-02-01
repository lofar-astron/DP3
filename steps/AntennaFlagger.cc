// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "AntennaFlagger.h"

#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "../base/FlagCounter.h"

namespace {
/**
 * Converts an input string with the format defined in preflagger's baseline
 * selection formats.
 * This can is further explained in
 * https://www.astron.nl/citt/DP3/steps/Description%20of%20baseline%20selection%20parameters.html
 *
 * Currently, the following selections are implemented;
 *   - Single antennas (e.g.  1,3,15,4)
 *   - Antenna ranges (e.g.  1~15)
 *   - Autocorrelations for a selection (e.g.  1~15&&& does all
 *       autocorrelations for the antennas in a selection)
 *  Multiple selections can be separated with a semicolon (;)
 * @param selection the input string to-be parsed.
 * @return a vector of pairs in the format (ant1, ant2)
 */
std::vector<std::pair<int, int>> ConvertSelection(
    const std::string& selection_string) {
  // The list of all baselines to be flagged
  std::vector<std::pair<int, int>> selection;

  // Split selection string on antenna-sets
  std::vector<std::string> antenna_sets;
  boost::split(antenna_sets, selection_string, boost::is_any_of(";"));

  for (std::string& antenna_set : antenna_sets) {
    // Split on antennas within the antenna-set.
    std::vector<std::string> antenna_strings;
    boost::split(antenna_strings, antenna_set, boost::is_any_of(","));

    // Check if correlation mode for this set of antennas.
    bool auto_correlation = false;
    bool cross_correlation = false;
    std::smatch matches;
    if (std::regex_search(antenna_set, matches, std::regex("&&&"))) {
      auto_correlation = true;
    } else if (std::regex_search(antenna_set, matches, std::regex("&&"))) {
      auto_correlation = true;
      cross_correlation = true;
    } else if (std::regex_search(antenna_set, matches, std::regex("&"))) {
      cross_correlation = true;
    } else {  // Nothing specified, same as &&
      auto_correlation = true;
      cross_correlation = true;
    }

    // Convert all antenna numbers from input string to a vector of ints.
    std::vector<int> antennas;
    for (const std::string& antenna : antenna_strings) {
      // Skip empty inputs
      if (antenna == "") {
        continue;
      }

      // This regex searches for a range, low~high (e.g. 0~20).
      std::regex range("(\\d+)~(\\d+)");
      std::smatch matches;

      try {
        if (std::regex_search(antenna, matches, range)) {
          // A range like 0~20 is in the string.
          const int high = std::stoi(matches[2]) + 1;  // High is inclusive
          for (int low = std::stoi(matches[1]); low < high; low++) {
            antennas.push_back(low);
          }
        } else {
          antennas.push_back(std::stoi(antenna));
        }
      } catch (std::invalid_argument& e) {
        std::cerr << "Could not parse: " << antenna << '\n';
        continue;
      }
    }

    if (auto_correlation) {
      selection.reserve(antennas.size());
      for (int antenna : antennas) {
        selection.emplace_back(antenna, antenna);
      }
    }

    if (cross_correlation) {
      for (int antenna1 : antennas) {
        selection.reserve(selection.size() + antenna1);
        for (int antenna2 = 0; antenna2 < antenna1; antenna2++) {
          selection.emplace_back(antenna1, antenna2);
        }
      }
    }
  }

  return selection;
}

xt::xtensor<bool, 2> Convert(const std::string& selection_string,
                             size_t n_antennas) {
  assert(n_antennas != 0);
  xt::xtensor<bool, 2> bl({n_antennas, n_antennas}, false);

  std::vector<std::pair<int, int>> selection =
      ConvertSelection(selection_string);

  for (std::pair<int, int> p : selection) {
    bl(p.first, p.second) = true;
    bl(p.second, p.first) = true;
  }

  return bl;
}
}  // anonymous namespace

namespace dp3 {
namespace steps {

void AntennaFlagger::finish() { getNextStep()->finish(); }

void AntennaFlagger::show(std::ostream& ostream) const {
  ostream << "AntennaFlagger " << name_;
  ostream << "\n  selection:   " << selection_string_;
}

AntennaFlagger::AntennaFlagger(const common::ParameterSet& parset,
                               const std::string& prefix)
    : name_(prefix),
      selection_string_(parset.getString(prefix + "selection", std::string())),
      do_detect_outliers_(parset.getBool(prefix + "detect_outliers", false)),
      antenna_flagging_sigma_(
          parset.getFloat(prefix + "antenna_flagging_sigma", 3)),
      antenna_flagging_maxiters_(
          parset.getInt(prefix + "antenna_flagging_maxiters", 5)),
      station_flagging_sigma_(
          parset.getFloat(prefix + "station_flagging_sigma", 2.5)),
      station_flagging_maxiters_(
          parset.getInt(prefix + "station_flagging_maxiters", 5)) {}

void AntennaFlagger::updateInfo(const base::DPInfo& info) {
  common::NSTimer::StartStop scoped_timer(initialization_timer_);

  Step::updateInfo(info);

  // Parse our selection string into a matrix of booleans.
  selection_ = Convert(selection_string_, info.nantenna());

  if (do_detect_outliers_) {
    std::string antenna_set(info.antennaSet());
    std::vector<std::string> antenna_names(info.antennaNames());
    size_t n_receivers_per_station = 0;
    size_t n_stations = 0;
    if (antenna_names[0].substr(0, 3) == "A12") {
      // AARTFAAC-12: all antennas in a station are used as individual receivers
      n_receivers_per_station = 48;
      const size_t n_antennas = info.nantenna();
      n_stations = n_antennas / n_receivers_per_station;
    } else {
      // e.g. LOFAR LBA or HBA: all antennas in a station are combined to form
      // one 'receiver'
      n_receivers_per_station = 1;
      n_stations = info.nantenna();
    }

    flagger_ = std::make_unique<dp3::antennaflagger::Flagger>(
        n_stations, n_receivers_per_station, info.nchan(), info.ncorr());
  }
}

bool AntennaFlagger::process(const base::DPBuffer& buffer) {
  computation_timer_.start();

  buffer_.copy(buffer);

  // Flags are in format: corr,chan,baseline.
  casacore::Cube<bool> flags = buffer_.GetCasacoreFlags();

  /*
   * Do additional flagging on the input data based on outlier statistics.
   */
  if (do_detect_outliers_) {
    computation_timer_.stop();
    flagging_timer_.start();
    const std::vector<std::complex<float>> data(
        buffer_.GetCasacoreData().begin(), buffer_.GetCasacoreData().end());
    flagger_->ComputeStats(data);

    // Find outliers both at antenna and station level
    const std::vector<size_t> flagged_antennas_1 = flagger_->FindBadAntennas(
        antenna_flagging_sigma_, antenna_flagging_maxiters_);
    const std::vector<size_t> flagged_antennas_2 = flagger_->FindBadStations(
        station_flagging_sigma_, station_flagging_maxiters_);

    // Combine vectors into a single set.
    std::set<size_t> flagged_antennas(flagged_antennas_2.begin(),
                                      flagged_antennas_2.end());
    flagged_antennas.insert(flagged_antennas_1.begin(),
                            flagged_antennas_1.end());

    // Flag all antennas
    for (int antenna : flagged_antennas) {
      for (size_t i = 0; i < info().nantenna(); i++) {
        selection_(antenna, i) = true;
        selection_(i, antenna) = true;
      }
    }

    flagging_timer_.stop();
    computation_timer_.start();
  }

  for (size_t correlation = 0; correlation < flags.nrow(); correlation++) {
    for (size_t channel = 0; channel < flags.ncolumn(); channel++) {
      for (size_t baseline = 0; baseline < flags.nplane(); baseline++) {
        const std::pair<size_t, size_t> antennas = dp3::common::ComputeBaseline(
            baseline, info().nantenna(), dp3::common::BaselineOrder::kRowMajor);
        flags(correlation, channel, baseline) =
            selection_(antennas.first, antennas.second);
      }
    }
  }

  buffer_.setFlags(flags);

  computation_timer_.stop();
  getNextStep()->process(buffer_);

  return true;
}

void AntennaFlagger::showTimings(std::ostream& ostream, double duration) const {
  double initializion_time = initialization_timer_.getElapsed();
  double flagging_time = flagging_timer_.getElapsed();
  double computation_time = computation_timer_.getElapsed();
  double total_time = initializion_time + flagging_time + computation_time;

  ostream << "  ";
  base::FlagCounter::showPerc1(ostream, total_time, duration);
  ostream << " AntennaFlagger " << name_ << "\n          ";
  base::FlagCounter::showPerc1(ostream, initializion_time, total_time);
  ostream << " of it spent in initialization.\n          ";
  base::FlagCounter::showPerc1(ostream, flagging_time, total_time);
  ostream << " of it spent in flagging.\n          ";
  base::FlagCounter::showPerc1(ostream, computation_time, total_time);
  ostream << " of it spent in setting flags.\n";
}

}  // namespace steps
}  // namespace dp3
