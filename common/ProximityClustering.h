// ProximityClustering.h: clustering to achieve speedup
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef PROXIMITY_CLUSTERING_H
#define PROXIMITY_CLUSTERING_H

#include <vector>
#include <utility>

namespace dp3 {
namespace common {

/**
 * Class that can cluster a source list into groups of sources
 * that are proximate to each other.
 */
class ProximityClustering {
 public:
  using NumType = double;  // can be changed to float if desired
  using Coordinate = std::pair<NumType, NumType>;

  ProximityClustering(const std::vector<Coordinate> &coordinates);
  std::vector<std::vector<std::size_t>> Group(NumType max_distance);

 private:
  Coordinate GetCoordinate(std::size_t i) const;
  NumType ClusterDistance(std::size_t i, std::size_t j) const;
  void GroupSource(std::size_t source_index, NumType max_distance);
  Coordinate Centroid(std::size_t i) const;
  static NumType EuclidDistance(Coordinate x1, Coordinate x2);

  std::vector<std::vector<std::size_t>> clusters_;
  const std::vector<Coordinate> &coordinates_;
};

}  // namespace common
}  // namespace dp3

#endif
