// ProximityClustering.cc: implementation of the proximity clustering class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ProximityClustering.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <cmath>
#include <cassert>

using std::size_t;

namespace dp3 {
namespace common {

ProximityClustering::ProximityClustering(
    const std::vector<ProximityClustering::Coordinate> &coordinates)
    : coordinates_(coordinates) {}

ProximityClustering::Coordinate ProximityClustering::GetCoordinate(
    size_t i) const {
  assert(i < coordinates_.size());
  return coordinates_[i];
}

ProximityClustering::NumType ProximityClustering::EuclidDistance(
    ProximityClustering::Coordinate x1, ProximityClustering::Coordinate x2) {
  return std::sqrt((x1.first - x2.first) * (x1.first - x2.first) +
                   (x1.second - x2.second) * (x1.second - x2.second));
}

ProximityClustering::Coordinate ProximityClustering::Centroid(size_t i) const {
  assert(!clusters_[i].empty());
  NumType x = 0;
  NumType y = 0;
  for (size_t j : clusters_[i]) {
    x += GetCoordinate(j).first;
    y += GetCoordinate(j).second;
  }
  x /= clusters_[i].size();
  y /= clusters_[i].size();
  return std::make_pair(x, y);
}

ProximityClustering::NumType ProximityClustering::ClusterDistance(
    size_t i, size_t j) const {
  Coordinate centroid1 = Centroid(i);
  Coordinate centroid2 = Centroid(j);
  return EuclidDistance(centroid1, centroid2);
}

void ProximityClustering::GroupSource(size_t source_index,
                                      NumType max_distance) {
  std::vector<std::vector<size_t>> new_clusters;
  for (size_t j = 0; j < source_index + 1; ++j) {
    new_clusters.push_back(clusters_[j]);
  }
  for (size_t j = source_index + 1; j < clusters_.size(); ++j) {
    if (ClusterDistance(source_index, j) < max_distance) {
      new_clusters[source_index].insert(new_clusters[source_index].end(),
                                        clusters_[j].begin(),
                                        clusters_[j].end());
    } else {
      new_clusters.push_back(clusters_[j]);
    }
  }
  clusters_ = std::move(new_clusters);
}

std::vector<std::vector<size_t>> ProximityClustering::Group(
    NumType max_distance) {
  clusters_.clear();
  for (size_t i = 0; i < coordinates_.size(); ++i) {
    clusters_.emplace_back(1, i);
  }

  size_t i = 0;
  while (i < clusters_.size()) {  // while loop because the number of clusters
                                  // decreases in an iteration
    GroupSource(i, max_distance);
    ++i;
  }
  return std::move(clusters_);
}

}  // namespace common
}  // namespace dp3
