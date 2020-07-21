// PrettyUnits.h - Print units in a human-readable way
//
// Copyright (C) 2008
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

#ifndef LOFAR_COMMON_PRETTY_UNITS_H
#define LOFAR_COMMON_PRETTY_UNITS_H

#include <string>

namespace DP3 {

/// \brief Print units in a human-readable way
class PrettyUnits : public std::string {
 protected:
  PrettyUnits(double value, const char *unit, unsigned precision);
};

/// \brief Print time in a human-readable way
class PrettyTime : public PrettyUnits {
 public:
  PrettyTime(double seconds = 0, unsigned precision = 3)
      : PrettyUnits(seconds, "s", precision) {}
};

/// \brief Print frequency in a human-readable way
class PrettyFrequency : public PrettyUnits {
 public:
  PrettyFrequency(double frequency = 0, unsigned precision = 3)
      : PrettyUnits(frequency, "Hz", precision) {}
};

}  // namespace DP3

#endif
