//# TypeNames.h: Return a string giving the type name to be stored in blobs
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
//# $Id: TypeNames.h 14057 2009-09-18 12:26:29Z diepen $

#ifndef LOFAR_COMMON_TYPENAMES_H
#define LOFAR_COMMON_TYPENAMES_H

// \file
// Return a string giving the type name to be stored in blobs

//# Includes
//#include "LofarTypes.h"

#include <complex>
#include <string>

namespace DP3
{
  // \addtogroup TypeNames
  //
  // These global functions return the name of the basic types.
  // They are meant to get the full id of a templated class when such an
  // object is stored in a blob.
  // As much as possible std::complex and builtin complex types get the same
  // name, so they can be read back from a blob in both ways.
  // <group>

  const std::string& typeName (const void*);
  const std::string& typeName (const bool*);
  const std::string& typeName (const char*);
  const std::string& typeName (const int8_t*);
  const std::string& typeName (const uint8_t*);
  const std::string& typeName (const int16_t*);
  const std::string& typeName (const uint16_t*);
  const std::string& typeName (const int32_t*);
  const std::string& typeName (const uint32_t*);
  const std::string& typeName (const int64_t*);
  const std::string& typeName (const uint64_t*);
  const std::string& typeName (const float*);
  const std::string& typeName (const double*);
  const std::string& typeName (const std::complex<float>*);
  const std::string& typeName (const std::complex<double>*);
#ifdef LOFAR_BUILTIN_COMPLEXINT
  const std::string& typeName (const std::complex<int16>*);
  const std::string& typeName (const std::complex<uint16>*);
#endif
  template<typename T> const std::string& typeName (T const* const*);

// </group>

}

// Include templated implementations.
#include "TypeNames.tcc"

#endif
