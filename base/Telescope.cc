// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Telescope.h"

#include <concepts>
#include <numeric>
#include <type_traits>

#include <EveryBeam/telescope/dish.h>
#include <EveryBeam/telescope/mwa.h>
#include <EveryBeam/telescope/phasedarray.h>

namespace dp3 {
namespace base {

std::vector<size_t> SelectStationIndices(
    const everybeam::telescope::Telescope& telescope,
    const std::vector<std::string>& station_names) {
  std::vector<size_t> station_to_msindex;
  if (telescope.IsHomogeneous()) {
    // All stations can be assumed to be identical
    station_to_msindex = {0};
    return station_to_msindex;
  }

  auto phased_array =
      dynamic_cast<const everybeam::telescope::PhasedArray*>(&telescope);
  assert(phased_array);

  // Copy only those stations for which the name matches.
  // Note: the order of the station names in both vectors match,
  // thus avoiding a nested loop.
  station_to_msindex.reserve(station_names.size());
  size_t station_idx = 0;
  for (size_t i = 0; i < telescope.GetNrStations(); ++i) {
    if (station_idx < station_names.size() &&
        phased_array->GetStation(i).GetName() == station_names[station_idx]) {
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
