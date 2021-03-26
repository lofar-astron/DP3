// PrettyUnits.cc - Print units in a human-readable way
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!

#include "PrettyUnits.h"

#include <sstream>
#include <iomanip>
#include <cmath>

namespace dp3 {
namespace common {

PrettyUnits::PrettyUnits(double value, const char *unit, unsigned precision) {
  static const char *prefixes = "yzafpnum kMGTPEZY";
  const char *prefix;

  if (value == 0.0)
    prefix = " ";
  else
    for (value *= 1e24, prefix = prefixes;
         fabs(value) >= 999.5 && prefix[1] != '\0'; prefix++)
      value /= 1000.0;

  std::stringstream stream;
  stream << std::setprecision(precision) << std::setw(precision + 1) << value;
  *static_cast<std::string *>(this) = stream.str() + ' ' + *prefix + unit;
}

}  // namespace common
}  // namespace dp3
