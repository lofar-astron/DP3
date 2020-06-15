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

#ifndef FITS_IO_CHECKER_H
#define FITS_IO_CHECKER_H

#include <string>

class FitsIOChecker
{
protected:
  static void checkStatus(int status, const std::string& filename);
  static void checkStatus(int status, const std::string& filename, const std::string& operation);
public:
  enum Unit {
    JanskyPerBeam,
    JanskyPerPixel,
    Jansky,
    Kelvin,
    MilliKelvin
  };
  static const char* UnitName(Unit unit) {
    switch(unit)
    {
      case JanskyPerBeam: return "Jansky/beam";
      case JanskyPerPixel: return "Jansky/pixel";
      case Jansky: return "Jansky";
      case Kelvin: return "Kelvin";
      case MilliKelvin: return "Milli-Kelvin";
    }
    return "";
  }
};

#endif
