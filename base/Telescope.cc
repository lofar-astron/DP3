// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Telescope.h"

#include <EveryBeam/telescope/phasedarray.h>

// Support both old and new return value types of PhasedArray::GetStation().
// TODO(AST-1001): Remove support for the old type in 2023.
namespace {
[[maybe_unused]] const std::string& GetStationName(
    const std::shared_ptr<const everybeam::Station>& station) {
  return station->GetName();
}

[[maybe_unused]] const std::string& GetStationName(
    const everybeam::Station& station) {
  return station.GetName();
}
}  // namespace

namespace dp3 {
namespace base {

std::vector<size_t> SelectStationIndices(
    const everybeam::telescope::Telescope& telescope,
    const std::vector<std::string>& station_names) {
  auto phased_array =
      dynamic_cast<const everybeam::telescope::PhasedArray*>(&telescope);
  if (phased_array == nullptr) {
    throw std::runtime_error(
        "Currently, only PhasedArray telescopes accepted as input, i.e. "
        "OSKAR or LOFAR. Support for other telescopes may become available "
        "soon.");
  }

  // Copy only those stations for which the name matches.
  // Note: the order of the station names in both vectors match,
  // thus avoiding a nested loop.
  std::vector<size_t> station_to_msindex;
  station_to_msindex.reserve(station_names.size());
  size_t station_idx = 0;
  for (size_t i = 0; i < phased_array->GetNrStations(); ++i) {
    if (station_idx < station_names.size() &&
        GetStationName(phased_array->GetStation(i)) ==
            station_names[station_idx]) {
      station_to_msindex.push_back(i);
      ++station_idx;
    }
  }

  if (station_idx != station_names.size()) {
    throw std::runtime_error(
        "SelectStationIndices: some stations miss the beam info");
  }
  return station_to_msindex;
}

}  // namespace base
}  // namespace dp3
