//# TypeNames.cc: Return a string giving the type name to be stored in blobs
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
//# $Id: TypeNames.cc 14057 2009-09-18 12:26:29Z diepen $

//# Always #include <lofar_config.h> first!

//# Includes
#include "TypeNames.h"

namespace DP3
{
  const std::string& typeName (const void*)
    {
      static std::string str ("unknown");
      return str;
    }

  const std::string& typeName (const bool*)
    {
      static std::string str ("bool");
      return str;
    }

  const std::string& typeName (const char*)
    {
      static std::string str ("char");
      return str;
    }

  const std::string& typeName (const int8_t*)
    {
      // This is also char (for backward compatibility).
      static std::string str ("char");
      return str;
    }

  const std::string& typeName (const uint8_t*)
    {
      static std::string str ("uchar");
      return str;
    }

  const std::string& typeName (const int16_t*)
    {
      static std::string str ("int16");
      return str;
    }

  const std::string& typeName (const uint16_t*)
    {
      static std::string str ("uint16");
      return str;
    }

  const std::string& typeName (const int32_t*)
    {
      static std::string str ("int32");
      return str;
    }

  const std::string& typeName (const uint32_t*)
    {
      static std::string str ("uint32");
      return str;
    }

  const std::string& typeName (const int64_t*)
    {
      static std::string str ("int64");
      return str;
    }

  const std::string& typeName (const uint64_t*)
    {
      static std::string str ("uint64");
      return str;
    }

  const std::string& typeName (const float*)
    {
      static std::string str ("float");
      return str;
    }

  const std::string& typeName (const double*)
    {
      static std::string str ("double");
      return str;
    }

  const std::string& typeName (const std::complex<float>*)
    {
      static std::string str ("fcomplex");
      return str;
    }
  const std::string& typeName (const std::complex<double>*)
    {
      static std::string str ("dcomplex");
      return str;
    }
}
