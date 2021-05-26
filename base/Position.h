// Position.h: A position on the celestial sphere.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_POSITION_H
#define DPPP_POSITION_H

#include <cstring>

namespace dp3 {
namespace base {

/// \brief A position on the celestial sphere.
/// @{

class Position {
 public:
  Position();

  /**
   * @brief Construct a new Position object
   *
   * @param alpha Hour angle in radians
   * @param delta Declination in radians
   */
  Position(double alpha, double delta);

  const double &operator[](size_t i) const;
  double &operator[](size_t i);

 private:
  double itsPosition[2];
};

/// @}

// -------------------------------------------------------------------------- //
// - Implementation: Position                                               - //
// -------------------------------------------------------------------------- //

inline const double &Position::operator[](size_t i) const {
  return itsPosition[i];
}

inline double &Position::operator[](size_t i) { return itsPosition[i]; }

}  // namespace base
}  // namespace dp3

#endif
