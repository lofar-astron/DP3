// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockTelescope.h"

namespace dp3::test {

MockPointResponse::MockPointResponse(const MockTelescope& telescope)
    : everybeam::pointresponse::PointResponse(&telescope, 0) {}

MockPointResponse::~MockPointResponse() {
  auto mock_telescope = static_cast<const MockTelescope&>(GetTelescope());
  BOOST_TEST(response_count_ == mock_telescope.ExpectedFrequencies().size(),
             "Unexpected number of PointResponse::Response calls.");
}

aocommon::MC2x2 MockPointResponse::Response(
    everybeam::BeamMode beam_mode, size_t station_idx, double freq,
    const everybeam::vector3r_t& direction, std::mutex*) {
  auto mock_telescope = static_cast<const MockTelescope&>(GetTelescope());
  BOOST_TEST((beam_mode == everybeam::BeamMode::kFull));
  BOOST_TEST(station_idx == 0);
  BOOST_TEST(freq == mock_telescope.ExpectedFrequencies()[response_count_]);
  const everybeam::vector3r_t& expected_direction =
      mock_telescope.ExpectedDirection();

  BOOST_CHECK_CLOSE(direction[0], expected_direction[0], 1.0e-6);
  BOOST_CHECK_CLOSE(direction[1], expected_direction[1], 1.0e-6);
  BOOST_CHECK_CLOSE(direction[2], expected_direction[2], 1.0e-6);

  ++response_count_;
  return aocommon::MC2x2(2.0, 0.0, 0.0, 2.0);
}

}  // namespace dp3::test
