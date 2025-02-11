// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>

#include <boost/test/unit_test.hpp>

#include "../../ProximityClustering.h"

using dp3::common::ProximityClustering;

namespace {
constexpr double kArcmin = M_PI / (180.0 * 60.0);
constexpr double kArcsec = M_PI / (180.0 * 60.0 * 60.0);
}  // namespace

BOOST_AUTO_TEST_SUITE(proximity_clustering)

BOOST_AUTO_TEST_CASE(small_set) {
  std::vector<std::pair<double, double>> coords{{10.0 * kArcmin, 1.0 * kArcmin},
                                                {0.0, 1.0 * kArcmin},
                                                {1.0 * kArcmin, 0.0},
                                                {10.0 * kArcmin, 0.0}};
  ProximityClustering pc(coords);
  std::vector<std::vector<size_t>> groups = pc.Group(2.0 * kArcmin);
  BOOST_REQUIRE_EQUAL(groups.size(), 2u);
  BOOST_REQUIRE_EQUAL(groups[0].size(), 2u);
  BOOST_REQUIRE_EQUAL(groups[1].size(), 2u);
  // Index 0 and 3 must be grouped together:
  BOOST_CHECK((groups[0][0] == 0u && groups[0][1] == 3u) ||
              (groups[0][0] == 3u && groups[0][1] == 0u) ||
              (groups[1][0] == 0u && groups[1][1] == 3u) ||
              (groups[1][0] == 3u && groups[1][1] == 0u));
  // Index 1 and 2 must be grouped together:
  BOOST_CHECK((groups[0][0] == 1u && groups[0][1] == 2u) ||
              (groups[0][0] == 2u && groups[0][1] == 1u) ||
              (groups[1][0] == 1u && groups[1][1] == 2u) ||
              (groups[1][0] == 2u && groups[1][1] == 1u));
}

BOOST_AUTO_TEST_CASE(large_set) {
  std::vector<std::pair<double, double>> coords;
  for (size_t i = 0; i != 100; ++i) {
    coords.emplace_back(0, 0);
    coords.emplace_back(1, 1);
    coords.emplace_back(2, 0);
  }
  ProximityClustering pc(coords);
  std::vector<std::vector<size_t>> groups = pc.Group(1.0);
  BOOST_REQUIRE_EQUAL(groups.size(), 3u);
  BOOST_REQUIRE_EQUAL(groups[0].size(), 100u);
  BOOST_REQUIRE_EQUAL(groups[1].size(), 100u);
  BOOST_REQUIRE_EQUAL(groups[2].size(), 100u);
  for (size_t i = 0; i != 100; ++i) {
    BOOST_CHECK_EQUAL(groups[0][i] % 3, 0u);
    BOOST_CHECK_EQUAL(groups[1][i] % 3, 1u);
    BOOST_CHECK_EQUAL(groups[2][i] % 3, 2u);
  }
}

BOOST_AUTO_TEST_CASE(around_ra_0) {
  // One source is placed just to before the RA 0 line and one
  // after, such that they should be clustered. The third one is
  // an unrelated distant third source.
  std::vector<std::pair<double, double>> coords{
      {2.0 * M_PI - kArcsec, 1.0}, {kArcsec, 1.0}, {100.0 * kArcsec, 1.0}};
  ProximityClustering pc(coords);
  std::vector<std::vector<size_t>> groups = pc.Group(10 * kArcsec);
  BOOST_REQUIRE_EQUAL(groups.size(), 2u);
  const size_t group_a = groups[0].size() == 2 ? 0 : 1;
  const size_t group_b = 1 - group_a;
  BOOST_REQUIRE_EQUAL(groups[group_a].size(), 2u);
  BOOST_REQUIRE_EQUAL(groups[group_b].size(), 1u);
  // Index 0 and 1 must be grouped together:
  BOOST_CHECK((groups[group_a][0] == 0u && groups[group_a][1] == 1u) ||
              (groups[group_a][0] == 1u && groups[group_a][1] == 0u));
}

BOOST_AUTO_TEST_SUITE_END()
