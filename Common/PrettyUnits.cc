//# PrettyUnits.cc - Print units in a human-readable way
//#
//# Copyright (C) 2008
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: PrettyUnits.cc 14057 2009-09-18 12:26:29Z diepen $

//# Always #include <lofar_config.h> first!

#include "PrettyUnits.h"

#include <sstream>
#include <iomanip>
#include <cmath>


namespace DP3 {

PrettyUnits::PrettyUnits(double value, const char *unit, unsigned precision)
{
  static const char *prefixes = "yzafpnum kMGTPEZY";
  const char	    *prefix;

  if (value == 0.0)
    prefix = " ";
  else
    for (value *= 1e24, prefix = prefixes; fabs(value) >= 999.5 && prefix[1] != '\0'; prefix ++)
      value /= 1000.0;

  std::stringstream stream;
  stream << std::setprecision(precision) << std::setw(precision + 1) << value;
  *static_cast<std::string *>(this) = stream.str() + ' ' + *prefix + unit;
}

}  // end namespace LOFAR
