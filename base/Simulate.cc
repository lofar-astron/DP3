// Simulate.cc: Simulate visibilities for a patch of sources.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "Simulate.h"

#include "../steps/PhaseShift.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/MatrixMath.h>

using casacore::Block;
using casacore::Matrix;
using casacore::Vector;

using dp3::steps::PhaseShift;

namespace dp3 {
namespace base {

std::vector<int> SetupUvwSplitting(unsigned int nant,
                                   const std::vector<int>& ant1,
                                   const std::vector<int>& ant2) {
  // Get the indices of the baselines needed to split the baseline UVWs into
  // station UVWs. They are in such an order that the UVW of a station is known
  // before used in another baseline to derive the UVW of the other station.
  // It can handle cases where baselines occur in disjoint station groups
  // like 0-1, 0-2, 1-2 and 3-4, 4-5, 5-6.
  // Note that the first station of a group gets UVW=0. All other station UVWs
  // are relative to it using the baseline UVWs.
  // Also note that nr of groups can be derived from the size of the returned
  // vector (because it contains no entry for the first antenna in a group).
  std::vector<int> uvw_bl;
  uvw_bl.reserve(nant);
  Block<bool> known(nant, false);
  unsigned int nset = 0;
  // Loop until all stations are set.
  while (nset < nant) {
    // Disjoint groups might exist, so keep a vector containing related antennae
    // which are members of the same group.
    std::vector<unsigned int> members(1);
    // Set first unset station as the reference station (which gets UVW=0).
    for (unsigned int i = 0; i < nant; ++i) {
      if (!known[i]) {
        members[0] = i;
        known[i] = true;
        ++nset;
        break;
      }
    }
    // Loop through all members in the group.
    // Note that new members can be appended in this loop.
    for (unsigned int j = 0; j < members.size(); ++j) {
      int refst = members[j];
      // Find all stations having a baseline with the reference station.
      for (unsigned int i = 0; i < ant1.size(); ++i) {
        int a1 = ant1[i];
        int a2 = ant2[i];
        // Only take baselines into account for which one station is known,
        // so the other can be derived from it.
        // The unknown station becomes member of the group.
        if (known[a1] != known[a2]) {
          if (a1 == refst) {
            uvw_bl.push_back(i);
            members.push_back(a2);
            known[a2] = true;
            ++nset;
          } else if (a2 == refst) {
            uvw_bl.push_back(-(i + 1));
            members.push_back(a1);
            known[a1] = true;
            ++nset;
          }
        }
      }
    }
  }
  return uvw_bl;
}

std::vector<size_t> GetBaselinesSortedByLength(
    const std::vector<int>& antennas1, const std::vector<int>& antennas2,
    const std::vector<std::array<double, 3>>& antenna_positions) {
  // sort baselines by length
  std::vector<double> baseline_length_squared(antennas1.size());
  for (size_t i = 0; i != baseline_length_squared.size(); ++i) {
    const int ant1 = antennas1[i];
    const int ant2 = antennas2[i];
    const double u = antenna_positions[ant2][0] - antenna_positions[ant1][0];
    const double v = antenna_positions[ant2][1] - antenna_positions[ant1][1];
    const double w = antenna_positions[ant2][2] - antenna_positions[ant1][2];
    baseline_length_squared[i] = u * u + v * v + w * w;
  }

  const auto compare_by_length = [&baseline_length_squared](size_t bl1,
                                                            size_t bl2) {
    return baseline_length_squared[bl1] > baseline_length_squared[bl2];
  };

  std::vector<size_t> bl_idx_sorted(antennas1.size());
  std::iota(bl_idx_sorted.begin(), bl_idx_sorted.end(), 0);
  std::sort(bl_idx_sorted.begin(), bl_idx_sorted.end(), compare_by_length);
  return bl_idx_sorted;
}

std::vector<size_t> GetBaselineSelection(
    size_t n_antennas, const std::vector<size_t>& sorted_baseline_ids,
    const std::vector<int>& antennas1, const std::vector<int>& antennas2) {
  // Initialize data structure to keep track of which baselines are selected,
  // to which group an antenna belongs and which antennas are in a group
  std::vector<size_t> bl_selection;
  constexpr int kUnset = -1;
  std::vector<int> ant_to_group_id_map(n_antennas, kUnset);
  // Each element is a group that contains antenna indices.
  std::vector<std::vector<int>> groups;

  for (size_t bl : sorted_baseline_ids) {
    const int ant1 = antennas1[bl];
    const int ant2 = antennas2[bl];
    const int antenna1_group = ant_to_group_id_map[ant1];
    const int antenna2_group = ant_to_group_id_map[ant2];
    if (antenna1_group == antenna2_group) {
      if (antenna1_group == kUnset) {
        // Both antennas are unkown
        // add both to a new group
        ant_to_group_id_map[ant1] = groups.size();
        ant_to_group_id_map[ant2] = groups.size();
        groups.push_back({ant1, ant2});
      } else {
        // Both antennas are known, and in the same group:
        // this baseline adds no new information
        continue;
      }
    } else {
      if (antenna1_group == kUnset) {
        // ant1 is unknown, add it to the group of ant2
        ant_to_group_id_map[ant1] = antenna2_group;
        groups[antenna2_group].push_back(ant1);
      } else if (antenna2_group == kUnset) {
        // ant2 is unknown, add it to the group of ant1
        ant_to_group_id_map[ant2] = antenna1_group;
        groups[antenna1_group].push_back(ant2);
      } else {
        // Both antennas are known, but in different groups
        // Add the group of ant2 to the group of ant1
        for (int ant : groups[antenna2_group]) {
          ant_to_group_id_map[ant] = antenna1_group;
        }
        groups[antenna1_group].insert(
            groups[antenna1_group].end(),
            groups[antenna2_group].begin(),  // OLD group of antenna2
            groups[antenna2_group].end());
        groups[antenna2_group].clear();
      }
    }
    // Add baseline to selection
    bl_selection.push_back(bl);
    // Exit early if selected baselines already form a single spanning tree
    if (bl_selection.size() == size_t(n_antennas - 1)) {
      break;
    }
  }
  return bl_selection;
}

std::vector<int> SetupUvwSplitting(
    unsigned int nant, const std::vector<int>& antennas1,
    const std::vector<int>& antennas2,
    const std::vector<std::array<double, 3>>& antenna_positions) {
  const std::vector<size_t> sorted_baselines =
      GetBaselinesSortedByLength(antennas1, antennas2, antenna_positions);

  const std::vector<size_t> bl_selection =
      GetBaselineSelection(nant, sorted_baselines, antennas1, antennas2);

  // Select the entries from vectors antenna1 and antenna2 that correspons to
  // selected baselines
  std::vector<int> antennas1_bl_selection;
  std::transform(bl_selection.begin(), bl_selection.end(),
                 std::back_inserter(antennas1_bl_selection),
                 [&antennas1](size_t bl) -> int { return antennas1[bl]; });
  std::vector<int> antennas2_bl_selection;
  std::transform(bl_selection.begin(), bl_selection.end(),
                 std::back_inserter(antennas2_bl_selection),
                 [&antennas2](size_t bl) -> int { return antennas2[bl]; });

  // Compute the indices for the selected baselines
  std::vector<int> bl_idx =
      SetupUvwSplitting(nant, antennas1_bl_selection, antennas2_bl_selection);

  // Map the indices in the selection to indices in the original baselines
  std::transform(bl_idx.begin(), bl_idx.end(), bl_idx.begin(),
                 [&bl_selection](int bl) -> int {
                   return bl >= 0 ? bl_selection[bl]
                                  : -bl_selection[-bl - 1] - 1;
                 });
  return bl_idx;
}

void SplitUvw(const std::vector<int>& baseline_indices,
              const std::vector<Baseline>& baselines,
              const DPBuffer::UvwType& uvw_bl,
              xt::xtensor<double, 2>& uvw_ant) {
  uvw_ant.fill(0.0);
  for (unsigned int i = 0; i < baseline_indices.size(); ++i) {
    int index = baseline_indices[i];
    if (index < 0) {
      // Ant2 is known.
      index = -index - 1;
      const size_t first_baseline = baselines[index].first;
      const size_t second_baseline = baselines[index].second;
      for (int j = 0; j < 3; ++j) {
        uvw_ant(first_baseline, j) =
            uvw_ant(second_baseline, j) - uvw_bl(index, j);
      }
    } else {
      // Ant1 is known.
      const size_t first_baseline = baselines[index].first;
      const size_t second_baseline = baselines[index].second;
      for (int j = 0; j < 3; ++j) {
        uvw_ant(second_baseline, j) =
            uvw_ant(first_baseline, j) + uvw_bl(index, j);
      }
    }
  }
}

void rotateUVW(const Direction& from, const Direction& to, size_t nUVW,
               double* uvw) {
  casacore::Matrix<double> oldUVW(3, 3);
  casacore::Matrix<double> newUVW(3, 3);
  PhaseShift::fillEulerMatrix(oldUVW, from);
  PhaseShift::fillEulerMatrix(newUVW, to);

  casacore::Matrix<double> tmp(
      casacore::product(casacore::transpose(newUVW), oldUVW));
  const double* R = tmp.data();

  for (size_t i = 0; i < 3 * nUVW; i += 3) {
    // Compute rotated UVW.
    double u = uvw[i + 0] * R[0] + uvw[i + 1] * R[3] + uvw[i + 2] * R[6];
    double v = uvw[i + 0] * R[1] + uvw[i + 1] * R[4] + uvw[i + 2] * R[7];
    double w = uvw[i + 0] * R[2] + uvw[i + 1] * R[5] + uvw[i + 2] * R[8];

    uvw[i + 0] = u;
    uvw[i + 1] = v;
    uvw[i + 2] = w;

    // Move to next station.
  }  // Stations.
}

}  // namespace base
}  // namespace dp3
