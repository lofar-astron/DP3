// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef COMMON_MEDIAN_H
#define COMMON_MEDIAN_H

#include <algorithm>
#include <vector>

namespace dp3 {
namespace common {

/** Find the median value in vector v
 * @tparm T Base type, typically float or double
 */
template <class T>
T Median(std::vector<T>& v) {
  if (v.empty()) {
    return static_cast<T>(0);
  }
  auto target = v.begin() + v.size() / 2;
  std::nth_element(v.begin(), target, v.end());
  if ((v.size() % 2) == 1) {
    return *target;
  } else {
    auto before = target - 1;
    std::nth_element(v.begin(), before, v.end());
    return (*target + *before) / static_cast<T>(2);
  }
}

}  // namespace common
}  // namespace dp3

#endif
