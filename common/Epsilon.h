// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Define comparison functions that use tolerances.
/// @author Maik Nijhuis

#ifndef COMMON_EPSILON_H
#define COMMON_EPSILON_H

#include <cmath>

namespace dp3 {
namespace common {

/** Check if two vectors are equal, given an epsilon.
 * @tparam T Base type, typically float or double.
 */
template <class T>
bool EpsilonEqual(const std::vector<T>& left, const std::vector<T>& right,
                  const T& epsilon) {
  if (left.size() != right.size()) {
    return false;
  }

  for (auto left_it = left.begin(), right_it = right.begin();
       left_it != left.end(); ++left_it, ++right_it) {
    if (std::abs(*left_it - *right_it) > epsilon) {
      return false;
    }
  }

  return true;
}

}  // namespace common
}  // namespace dp3

#endif
