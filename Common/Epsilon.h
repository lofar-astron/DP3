// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Define comparison functions that use tolerances.
/// @author Maik Nijhuis

#ifndef COMMON_EPSILON_H
#define COMMON_EPSILON_H

#include <cmath>

namespace DP3 {

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

}  // namespace DP3

#endif
