// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../ProximityClustering.h"

BOOST_AUTO_TEST_SUITE(proximity_clustering)

BOOST_AUTO_TEST_CASE(small_set) {
  using DP3::ProximityClustering;

  std::vector<std::pair<double, double>> coords{
      {10.0, 1.0}, {0.0, 1.0}, {1.0, 0.0}, {10.0, 0.0}};
  ProximityClustering pc(coords);
  std::vector<std::vector<size_t>> groups = pc.Group(2.0);
  BOOST_CHECK_EQUAL(groups.size(), 2);
  BOOST_CHECK_EQUAL(groups[0].size(), 2);
  BOOST_CHECK_EQUAL(groups[1].size(), 2);
  // Index 0 and 3 must be grouped together:
  BOOST_CHECK((groups[0][0] == 0 && groups[0][1] == 3) ||
              (groups[0][0] == 3 && groups[0][1] == 0) ||
              (groups[1][0] == 0 && groups[1][1] == 3) ||
              (groups[1][0] == 3 && groups[1][1] == 0));
  // Index 1 and 2 must be grouped together:
  BOOST_CHECK((groups[0][0] == 1 && groups[0][1] == 2) ||
              (groups[0][0] == 2 && groups[0][1] == 1) ||
              (groups[1][0] == 1 && groups[1][1] == 2) ||
              (groups[1][0] == 2 && groups[1][1] == 1));
}

BOOST_AUTO_TEST_CASE(large_set) {
  using DP3::ProximityClustering;

  std::vector<std::pair<double, double>> coords;
  for (size_t i = 0; i != 100; ++i) {
    coords.emplace_back(0, 0);
    coords.emplace_back(1, 1);
    coords.emplace_back(2, 0);
  }
  ProximityClustering pc(coords);
  std::vector<std::vector<size_t>> groups = pc.Group(1.0);
  BOOST_CHECK_EQUAL(groups.size(), 3);
  BOOST_CHECK_EQUAL(groups[0].size(), 100);
  BOOST_CHECK_EQUAL(groups[1].size(), 100);
  BOOST_CHECK_EQUAL(groups[2].size(), 100);
  for (size_t i = 0; i != 100; ++i) {
    BOOST_CHECK_EQUAL(groups[0][i] % 3, 0);
    BOOST_CHECK_EQUAL(groups[1][i] % 3, 1);
    BOOST_CHECK_EQUAL(groups[2][i] % 3, 2);
  }
}

BOOST_AUTO_TEST_SUITE_END()
