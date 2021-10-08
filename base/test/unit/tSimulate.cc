// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Unit tests for functions in Simulate.cc.
/// @author Sebastiaan van der Tol

#include <boost/test/unit_test.hpp>
#include <algorithm>

#include "../../Simulate.h"

using dp3::base::nsetupSplitUVW;

BOOST_AUTO_TEST_SUITE(simulate)

BOOST_AUTO_TEST_CASE(nsetupsplituvw) {
  //
  // Simple test for nsetupSplitUVW with antenna position info
  //
  // The antenna layout is as follows
  //
  // 1  3  5---6
  // |  |  |  /
  // |  |  | /
  // |  |  |/
  // 0  2--4
  //
  // One of the baselines in the 4-5-6 cycle is redundant
  // Because baseline 5-6 is the shortest in the cycle, it should not be
  // selected All other baselines are needed to make the spanning trees
  //

  std::vector<std::array<double, 3>> antenna_pos = {
      {0, 0, 0}, {0, 4, 0}, {2, 0, 0}, {2, 4, 0},
      {4, 0, 0}, {4, 4, 0}, {7, 4, 0},
  };
  std::vector<int> antenna1 = {0, 3, 2, 4, 5, 4};
  std::vector<int> antenna2 = {1, 2, 4, 5, 6, 6};

  size_t nant = antenna_pos.size();

  std::vector<int> bl_indices =
      nsetupSplitUVW(nant, antenna1, antenna2, antenna_pos);

  // Baseline nr. 4 (antennas 5-6) should be excluded
  std::vector<int> good_bl_indices = {0, 1, 2, 3, 5};

  BOOST_CHECK_EQUAL(bl_indices.size(), good_bl_indices.size());

  // All good_indices must be present in the bl_indices found by nsetupSplitUVW
  for (int good_bl_idx : good_bl_indices) {
    // Reversed baselines are represented by negative indices
    // These are a matches too
    auto indices_match = [good_bl_idx](int idx) -> bool {
      return ((good_bl_idx == idx) || (good_bl_idx == -(idx + 1)));
    };
    BOOST_CHECK(std::find_if(bl_indices.begin(), bl_indices.end(),
                             indices_match) != bl_indices.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()
