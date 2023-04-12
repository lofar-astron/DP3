// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <dp3/common/Fields.h>

#include <string>
#include <ostream>
#include <vector>

namespace dp3 {
namespace common {

std::ostream& operator<<(std::ostream& stream, const Fields& fields) {
  std::vector<std::string_view> strings;
  if (fields.Data()) strings.push_back("data");
  if (fields.Flags()) strings.push_back("flags");
  if (fields.Weights()) strings.push_back("weights");
  if (fields.Uvw()) strings.push_back("uvw");

  stream << "[";
  for (size_t i = 0; i < strings.size(); ++i) {
    if (i != 0) stream << ", ";
    stream << strings[i];
  }
  stream << "]";

  return stream;
}

}  // namespace common
}  // namespace dp3
