// Position.h: A position on the celestial sphere.
//
// Copyright (C) 2012
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

#ifndef DPPP_POSITION_H
#define DPPP_POSITION_H

#include <cstring>

namespace DP3 {
namespace DPPP {

/// \brief A position on the celestial sphere.
/// @{

class Position {
 public:
  Position();
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

}  // namespace DPPP
}  // namespace DP3

#endif
