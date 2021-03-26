// TypeNames.tcc: Return a string giving the type name to be stored in blobs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later



#ifndef COMMON_TYPENAMES_TCC
#define COMMON_TYPENAMES_TCC

#include "TypeNames.h"

namespace dp3 {
namespace common {
  
template<typename T>
const std::string& typeName (T const* const*)
{
  static std::string str ("array<" + typeName((const T*)0) + ">");
  return str;
}

}
}

#endif
