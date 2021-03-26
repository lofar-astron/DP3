// ParameterRecord.cc: A record of parameter values
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParameterRecord.h"

#include <cstdio>
#include <ostream>
#include <string>

namespace dp3 {
namespace common {

std::ostream& operator<<(std::ostream& os, const ParameterRecord& prec) {
  bool first = true;
  os << '{';
  for (ParameterRecord::const_iterator iter = prec.begin(); iter != prec.end();
       ++iter) {
    if (first) {
      first = false;
    } else {
      os << ", ";
    }
    os << '\'' << iter->first << "': " << iter->second;
  }
  os << '}';
  return os;
}

bool ParameterRecord::getRecursive(const std::string& key,
                                   ParameterValue& value) const {
  const_iterator iter = find(key);
  if (iter != end()) {
    value = iter->second;
    return true;
  }
  // Try to find the key in possible ParmRecords.
  // Strip the last part from the key.
  std::string::size_type pos = key.rfind('.');
  while (pos != std::string::npos) {
    const_iterator iter = find(key.substr(0, pos));
    if (iter != end()) {
      ParameterValue pv(iter->second);
      if (pv.isRecord() &&
          pv.getRecord().getRecursive(key.substr(pos + 1), value)) {
        return true;
      }
    }
    if (pos == 0) return false;
    pos = key.rfind('.', pos - 1);
  }
  return false;
}

}  // namespace common
}  // namespace dp3
