// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Telescope.h"

#include <numeric>

#include <EveryBeam/telescope/dish.h>
#include <EveryBeam/telescope/mwa.h>
#include <EveryBeam/telescope/phasedarray.h>

namespace dp3 {
namespace base {

bool IsHomogeneous(const everybeam::telescope::Telescope& telescope) {
  // It would be better if EveryBeam provides this information ; a dish isn't
  // necessarily always homogenous (in future implementations).
  return dynamic_cast<const everybeam::telescope::Dish*>(&telescope) ||
         dynamic_cast<const everybeam::telescope::MWA*>(&telescope);
}

std::vector<size_t> SelectStationIndices(
    const everybeam::telescope::Telescope& telescope,
    const std::vector<std::string>& station_names) {
  std::vector<size_t> station_to_msindex;
  if (IsHomogeneous(telescope)) {
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
