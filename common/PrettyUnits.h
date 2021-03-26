// PrettyUnits.h - Print units in a human-readable way
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PRETTY_UNITS_H
#define LOFAR_COMMON_PRETTY_UNITS_H

#include <string>

namespace dp3 {
namespace common {

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

}  // namespace common
}  // namespace dp3

#endif
