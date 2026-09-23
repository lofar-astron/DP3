// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_TELESCOPE_H_
#define DP3_BASE_TELESCOPE_H_

#include <EveryBeam/everybeam.h>
#include <EveryBeam/telescope.h>

namespace dp3 {
namespace base {

/**
 * Retrieve the everybeam telescope from a Measurement Set.
 */
inline std::unique_ptr<everybeam::Telescope> GetTelescope(
    const casacore::MeasurementSet& ms,
    const everybeam::ElementResponseModel element_response_model,
    bool use_channel_frequency, const std::string& coefficients_file) {
  everybeam::Options options;
  options.element_response_model = element_response_model;
  options.use_channel_frequency = use_channel_frequency;
  options.coeff_path = coefficients_file;
  std::unique_ptr<everybeam::Telescope> telescope =
      everybeam::LoadTelescope(ms, options);
  return telescope;
}

/**
 * Find stations in a telescope by name and return their indices.
 * @param telescope The telescope, which contains antennae / stations.
 * @param station_names A list of station names. The order of the names in this
 *        list should match the order in which they occur in the telescope.
 * @return The indices corresponding to the station names. Because of the
 *         ordering restriction, the list always has increasing indices only.
 */
std::vector<size_t> GetStationIndices(
    const everybeam::Telescope& telescope,
    const std::vector<std::string>& station_names);

}  // namespace base
}  // namespace dp3

#endif
