//# TypeNames.tcc: Return a string giving the type name to be stored in blobs
//#
//# Copyright (C) 2003
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
//# $Id: TypeNames.tcc 14057 2009-09-18 12:26:29Z diepen $


#ifndef COMMON_TYPENAMES_TCC
#define COMMON_TYPENAMES_TCC

//# Includes
#include "TypeNames.h"

namespace DP3
{
  template<typename T>
  const std::string& typeName (T const* const*)
    {
      static std::string str ("array<" + typeName((const T*)0) + ">");
      return str;
    }
}

#endif
