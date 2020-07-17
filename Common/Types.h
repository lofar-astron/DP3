// Types.h: Define common types.
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
/// @brief Define common types.
/// @author Maik Nijhuis

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <casacore/casa/version.h>

namespace DP3 {

#if CASACORE_MAJOR_VERSION<3 || (CASACORE_MAJOR_VERSION==3 && CASACORE_MINOR_VERSION<4)
  typedef unsigned int rownr_t;
#else
  typedef casacore::rownr_t rownr_t;
#endif

}

#endif