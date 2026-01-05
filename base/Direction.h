// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_DIRECTION_H_
#define DP3_BASE_DIRECTION_H_

#include <cstring>

namespace dp3 {
namespace base {

/// \brief A direction on the celestial sphere.
/// @{

struct Direction {
  constexpr Direction() : ra(0.0), dec(0.0) {}

  /**
   * @brief Construct a new Direction object
   *
   * @param ra Right ascension in radians
   * @param dec Declination in radians
   */
  constexpr Direction(double _ra, double _dec) : ra(_ra), dec(_dec) {}

  double ra;   ///< Right ascension in radians
  double dec;  ///< Declination in radians
};

/// @}

}  // namespace base
}  // namespace dp3

#endif  // DP3_BASE_DIRECTION_H_
