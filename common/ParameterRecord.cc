// ParameterRecord.cc: A record of parameter values
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParameterRecord.h"

#include <ostream>

namespace dp3 {
namespace common {

std::ostream& operator<<(std::ostream& os, const ParameterRecord& record) {
  bool first = true;
  os << '{';
  for (const std::pair<std::string, ParameterValue>& entry : record) {
    if (first) {
      first = false;
    } else {
      os << ", ";
    }
    os << '\'' << entry.first << "': " << entry.second;
  }
  os << '}';
  return os;
}

}  // namespace common
}  // namespace dp3
