// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../UVWCalculator.h"

#include <boost/test/unit_test.hpp>

#include <vector>

using dp3::base::UVWCalculator;

namespace {
std::array<double, 3> GetUvw(UVWCalculator& calc, unsigned int ant1,
                             unsigned int ant2) {
  const double kTime = 44362.42;
  const std::array<double, 3> forward = calc.getUVW(ant1, ant2, kTime);
  const std::array<double, 3> backward = calc.getUVW(ant2, ant1, kTime);
  BOOST_CHECK_CLOSE(forward[0], -backward[0], 1.0e-6);
  BOOST_CHECK_CLOSE(forward[1], -backward[1], 1.0e-6);
  BOOST_CHECK_CLOSE(forward[2], -backward[2], 1.0e-6);
  return forward;
}

/**
 * Test that UVWCalculator::getUVW() satisfies:
 * - getUVW(x, y) == -getUVW(y, x)
 * - getUVW(x, y) + getUVW(y, z) == getUVW(x, z)
 */
void TestRelativeUvw(const casacore::MDirection& phase_direction) {
  const std::vector<double> array_values{0, 0, 1};
  const std::vector<std::vector<double>> station_values{
      {0, 0, 1}, {1, 0, 0}, {0, 2, 3}};

  const casacore::MPosition array_position(
      casacore::Quantum<casacore::Vector<double>>(array_values, "m"),
      casacore::MPosition::ITRF);
  std::vector<casacore::MPosition> station_positions;
  for (const std::vector<double>& values : station_values) {
    station_positions.emplace_back(
        casacore::Quantum<casacore::Vector<double>>(values, "m"),
        casacore::MPosition::ITRF);
  }

  UVWCalculator calc(phase_direction, array_position, station_positions);
  const std::array<double, 3> uvw_0_1 = GetUvw(calc, 0, 1);
  const std::array<double, 3> uvw_1_2 = GetUvw(calc, 1, 2);
  const std::array<double, 3> uvw_0_2 = GetUvw(calc, 0, 2);

  BOOST_CHECK_CLOSE(uvw_0_1[0] + uvw_1_2[0], uvw_0_2[0], 1.0e-6);
  BOOST_CHECK_CLOSE(uvw_0_1[1] + uvw_1_2[1], uvw_0_2[1], 1.0e-6);
  BOOST_CHECK_CLOSE(uvw_0_1[2] + uvw_1_2[2], uvw_0_2[2], 1.0e-6);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(uvwcalculator)

BOOST_AUTO_TEST_CASE(relative_uvw) {
  const casacore::MDirection phase_direction(casacore::Quantity(0, "deg"),
                                             casacore::Quantity(90, "deg"),
                                             casacore::MDirection::J2000);
  TestRelativeUvw(phase_direction);
}

BOOST_AUTO_TEST_CASE(relative_uvw_moving_phasedirection) {
  TestRelativeUvw(casacore::MDirection::makeMDirection("Sun"));
}

BOOST_AUTO_TEST_SUITE_END()