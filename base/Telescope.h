// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_TELESCOPE_H_
#define DP3_BASE_TELESCOPE_H_

#include <EveryBeam/load.h>
#include <EveryBeam/telescope/phasedarray.h>
#include <EveryBeam/telescope/dish.h>

namespace dp3 {
namespace base {

/**
 * Retrieve the everybeam telescope from a Measurement Set.
 */
inline std::unique_ptr<everybeam::telescope::Telescope> GetTelescope(
    const std::string& ms_name,
    const everybeam::ElementResponseModel element_response_model,
    bool use_channel_frequency) {
  everybeam::Options options;
  options.element_response_model = element_response_model;
  options.use_channel_frequency = use_channel_frequency;
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      everybeam::Load(ms_name, options);
  return telescope;
}

inline bool IsPhasedArray(const everybeam::telescope::Telescope& telescope) {
  auto phased_array =
      dynamic_cast<const everybeam::telescope::PhasedArray*>(&telescope);
  return phased_array != nullptr;
}

inline bool IsDish(const everybeam::telescope::Telescope& telescope) {
  auto dish = dynamic_cast<const everybeam::telescope::Dish*>(&telescope);
  return dish != nullptr;
}

/**
 * Find stations in a telescope by name and return their indices.
 * @param telescope The telescope, which contains antennae / stations.
 * @param station_names A list of station names. The order of the names in this
 *        list should match the order in which they occur in the telescope.
 * @return The indices corresponding to the station names. Because of the
 *         ordering restriction, the list always has increasing indices only.
 */
std::vector<size_t> SelectStationIndices(
    const everybeam::telescope::Telescope& telescope,
    const std::vector<std::string>& station_names);

}  // namespace base
}  // namespace dp3

#endif
