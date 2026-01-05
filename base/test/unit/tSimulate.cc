// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Unit tests for functions in Simulate.cc.
/// @author Sebastiaan van der Tol

#include "base/Simulate.h"

#include <algorithm>
#include <array>
#include <random>
#include <set>

#include <boost/test/unit_test.hpp>

using dp3::base::GetBaselineSelection;
using dp3::base::GetBaselinesSortedByLength;
using dp3::base::SetupUvwSplitting;
using dp3::base::SplitUvw;

namespace {

std::vector<std::array<double, 3>> GenerateRandomAntennaPositions(
    size_t n_antennas) {
  std::mt19937 generator;
  std::uniform_real_distribution<double> position_distribution(0.0, 1000.0);
  std::vector<std::array<double, 3>> antenna_positions(n_antennas);
  for (std::array<double, 3>& position : antenna_positions) {
    for (double& element : position) element = position_distribution(generator);
  }
  // Antenna 0 will be used as the reference antenna, so make it zero
  // for ease.
  antenna_positions[0] = {0.0, 0.0, 0.0};
  return antenna_positions;
}

std::vector<std::array<double, 3>> GenerateIncreasingAntennaPositions(
    size_t n_antennas) {
  std::vector<std::array<double, 3>> antenna_positions(n_antennas);
  for (size_t i = 0; i != antenna_positions.size(); ++i) {
    const double x = i;
    antenna_positions[i] = {x, x, x};
  }
  return antenna_positions;
}

void TestMwaCase(bool randomize_antenna_positions) {
  constexpr size_t kNAntennas = 128;
  constexpr size_t kNBaselines = kNAntennas * (kNAntennas + 1) / 2;
  std::vector<int> antenna1;
  std::vector<int> antenna2;
  std::vector<std::pair<size_t, size_t>> baselines;
  for (size_t a1 = 0; a1 != kNAntennas; ++a1) {
    for (size_t a2 = a1; a2 != kNAntennas; ++a2) {
      antenna1.emplace_back(a1);
      antenna2.emplace_back(a2);
      baselines.emplace_back(a1, a2);
    }
  }

  const std::vector<std::array<double, 3>> antenna_positions =
      randomize_antenna_positions
          ? GenerateRandomAntennaPositions(kNAntennas)
          : GenerateIncreasingAntennaPositions(kNAntennas);

  const std::vector<size_t> sorted_baselines =
      GetBaselinesSortedByLength(antenna1, antenna2, antenna_positions);
  const std::vector<size_t> baseline_selection =
      GetBaselineSelection(kNAntennas, sorted_baselines, antenna1, antenna2);

  // Every antenna should be reachable from the first baseline
  std::set<int> reached;
  reached.insert(antenna1[baseline_selection[0]]);
  reached.insert(antenna2[baseline_selection[0]]);
  bool repeat_search;
  do {
    repeat_search = false;
    for (int antenna : reached) {
      for (size_t baseline : baseline_selection) {
        if (antenna == antenna1[baseline] &&
            reached.count(antenna2[baseline]) == 0) {
          reached.insert(antenna2[baseline]);
          repeat_search = true;
          break;
        }
        if (antenna == antenna2[baseline] &&
            reached.count(antenna1[baseline]) == 0) {
          reached.insert(antenna1[baseline]);
          repeat_search = true;
          break;
        }
      }
      if (repeat_search)
        break;  // Exit the loop over 'reached', since it changed.
    }
  } while (repeat_search);
  BOOST_TEST(reached.size() == kNAntennas);

  // We use the ITRF positions of the antenna also as their UVW, which of course
  // is not realistic but should work and makes checking a bit more easy.
  dp3::base::DPBuffer::UvwType baseline_uvws({kNBaselines, 3});
  for (size_t baseline = 0; baseline != kNBaselines; ++baseline) {
    const std::array<double, 3>& position1 =
        antenna_positions[antenna1[baseline]];
    const std::array<double, 3>& position2 =
        antenna_positions[antenna2[baseline]];
    baseline_uvws(baseline, 0) = position2[0] - position1[0];
    baseline_uvws(baseline, 1) = position2[1] - position1[1];
    baseline_uvws(baseline, 2) = position2[2] - position1[2];
  }

  const std::vector<int> indices =
      SetupUvwSplitting(kNAntennas, antenna1, antenna2, antenna_positions);

  xt::xtensor<double, 2> antenna_uvws({kNAntennas, 3});
  SplitUvw(indices, baselines, baseline_uvws, antenna_uvws);

  // Check if calculated UVWs match the positions
  for (size_t antenna = 0; antenna != kNAntennas; ++antenna) {
    const std::array<double, 3>& reference = antenna_positions[antenna];
    BOOST_CHECK_CLOSE_FRACTION(reference[0], antenna_uvws(antenna, 0), 1.0e-6);
    BOOST_CHECK_CLOSE_FRACTION(reference[1], antenna_uvws(antenna, 1), 1.0e-6);
    BOOST_CHECK_CLOSE_FRACTION(reference[2], antenna_uvws(antenna, 2), 1.0e-6);
  }
}

}  // namespace

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
      SetupUvwSplitting(nant, antenna1, antenna2, antenna_pos);

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

// BOOST_DATA_TEST_CASE is not used here because it names test cases like
// _0 and _1, which makes it harder to analyse and debug.
BOOST_AUTO_TEST_CASE(uvws_mwa_case_simple) { TestMwaCase(false); }

BOOST_AUTO_TEST_CASE(uvws_mwa_case_randomized_positions) { TestMwaCase(true); }

BOOST_AUTO_TEST_SUITE_END()
