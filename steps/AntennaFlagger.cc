// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "AntennaFlagger.h"

#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <aocommon/xt/span.h>
#include <boost/algorithm/string.hpp>
#include <xtensor/xview.hpp>

#include "../base/FlagCounter.h"

namespace {
bool HasBaselineOrder(const dp3::base::DPInfo& info,
                      dp3::common::BaselineOrder baseline_order) {
  const std::vector<int>& ant1 = info.getAnt1();
  const std::vector<int>& ant2 = info.getAnt2();

  for (size_t bl = 0; bl < info.nbaselines(); ++bl) {
    const std::pair<int, int> antenna_pair =
        dp3::common::ComputeBaseline(bl, info.nantenna(), baseline_order);
    if (antenna_pair.first != ant1[bl] || antenna_pair.second != ant2[bl]) {
      return false;
    }
  }

  return true;
}
}  // namespace

namespace dp3 {
namespace steps {

void AntennaFlagger::finish() { getNextStep()->finish(); }

void AntennaFlagger::show(std::ostream& ostream) const {
  ostream << "AntennaFlagger " << name_;
}

AntennaFlagger::AntennaFlagger(const common::ParameterSet& parset,
                               const std::string& prefix)
    : name_(prefix),
      antenna_flagging_sigma_(
          parset.getFloat(prefix + "antenna_flagging_sigma", 3)),
      antenna_flagging_max_iterations_(
          parset.getInt(prefix + "antenna_flagging_max_iterations", 5)),
      station_flagging_sigma_(
          parset.getFloat(prefix + "station_flagging_sigma", 2.5)),
      station_flagging_max_iterations_(
          parset.getInt(prefix + "station_flagging_max_iterations", 5)) {}

void AntennaFlagger::updateInfo(const base::DPInfo& info) {
  common::NSTimer::StartStop scoped_timer(initialization_timer_);

  Step::updateInfo(info);

  if (dp3::common::ComputeNBaselines(info.nantenna()) != info.nbaselines()) {
    throw std::runtime_error(
        "AntennaFlagger requires a baseline for all antenna pairs.");
  }

  if (HasBaselineOrder(info, common::BaselineOrder::kColumnMajor)) {
    baseline_order_ = common::BaselineOrder::kColumnMajor;
  } else if (HasBaselineOrder(info, common::BaselineOrder::kRowMajor)) {
    baseline_order_ = common::BaselineOrder::kRowMajor;
  } else {
    throw std::runtime_error("Unsupported baseline order.");
  }

  size_t n_antennas_per_station = 0;
  size_t n_stations = 0;
  const std::string kAartfaac12Prefix = "A12";
  if (info.antennaNames()[0].substr(0, 3) == kAartfaac12Prefix) {
    // AARTFAAC-12: all antennas in a station used individually
    n_antennas_per_station = 48;
    if (info.nantenna() % n_antennas_per_station) {
      throw std::runtime_error(
          "Number of antennas not divisible by 48 as expected for "
          "AARTFAAC-12.");
    }
    const size_t n_antennas = info.nantenna();
    n_stations = n_antennas / n_antennas_per_station;
  } else {
    // e.g. LOFAR LBA or HBA: all antennas in a station are combined
    n_antennas_per_station = 1;
    n_stations = info.nantenna();
  }

  flagger_ = std::make_unique<dp3::antennaflagger::Flagger>(
      n_stations, n_antennas_per_station, info.nchan(), info.ncorr());
}

bool AntennaFlagger::process(std::unique_ptr<base::DPBuffer> buffer) {
  computation_timer_.start();
  flagger_->ComputeStats(buffer->GetData(), baseline_order_);

  // Find outlier antennas
  const xt::xtensor<int, 1> antenna_flags =
      flagger_->FindBadAntennas(antenna_flagging_sigma_,
                                antenna_flagging_max_iterations_) |
      flagger_->FindBadStations(station_flagging_sigma_,
                                station_flagging_max_iterations_);
  const std::vector<size_t> flagged_antennas = xt::where(antenna_flags)[0];
  computation_timer_.stop();

  // Set baseline flags
  flagging_timer_.start();

  if (buffer->GetFlags().size() == 0) {
    buffer->GetFlags().resize({getInfoOut().nbaselines(), getInfoOut().nchan(),
                               getInfoOut().ncorr()});
    buffer->GetFlags().fill(false);
  }
  base::DPBuffer::FlagsType& flags = buffer->GetFlags();

  for (size_t antenna1 : flagged_antennas) {
    for (size_t antenna2 = 0; antenna2 < getInfoOut().nantenna(); ++antenna2) {
      const size_t bl = common::ComputeBaselineIndex(
          antenna1, antenna2, getInfoOut().nantenna(), baseline_order_);
      xt::view(flags, bl, xt::all(), xt::all()) = true;
    }
  }

  flagging_timer_.stop();

  getNextStep()->process(std::move(buffer));

  return true;
}

void AntennaFlagger::showTimings(std::ostream& ostream, double duration) const {
  double initializion_time = initialization_timer_.getElapsed();
  double computation_time = computation_timer_.getElapsed();
  double flagging_time = flagging_timer_.getElapsed();
  double total_time = initializion_time + computation_time + flagging_time;

  ostream << "  ";
  base::FlagCounter::showPerc1(ostream, total_time, duration);
  ostream << " AntennaFlagger " << name_ << "\n          ";
  base::FlagCounter::showPerc1(ostream, initializion_time, total_time);
  ostream << " of it spent in initialization.\n          ";
  base::FlagCounter::showPerc1(ostream, computation_time, total_time);
  ostream << " of it spent in computing statistics and flags.\n          ";
  base::FlagCounter::showPerc1(ostream, flagging_time, total_time);
  ostream << " of it spent in setting flags.\n";
}

}  // namespace steps
}  // namespace dp3
