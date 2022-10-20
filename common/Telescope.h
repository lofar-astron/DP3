// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_TELESCOPE_H_
#define DP3_COMMON_TELESCOPE_H_

#include <EveryBeam/load.h>

namespace dp3 {
namespace common {

/// Retrieve the everybeam telescope from a Measurement Set.
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

}  // namespace common
}  // namespace dp3

#endif
