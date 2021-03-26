// TypeNames.cc: Return a string giving the type name to be stored in blobs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!

#include "TypeNames.h"

namespace dp3 {
namespace common {

const std::string& typeName(const void*) {
  static std::string str("unknown");
  return str;
}

const std::string& typeName(const bool*) {
  static std::string str("bool");
  return str;
}

const std::string& typeName(const char*) {
  static std::string str("char");
  return str;
}

const std::string& typeName(const int8_t*) {
  // This is also char (for backward compatibility).
  static std::string str("char");
  return str;
}

const std::string& typeName(const uint8_t*) {
  static std::string str("uchar");
  return str;
}

const std::string& typeName(const int16_t*) {
  static std::string str("int16");
  return str;
}

const std::string& typeName(const uint16_t*) {
  static std::string str("uint16");
  return str;
}

const std::string& typeName(const int32_t*) {
  static std::string str("int32");
  return str;
}

const std::string& typeName(const uint32_t*) {
  static std::string str("uint32");
  return str;
}

const std::string& typeName(const int64_t*) {
  static std::string str("int64");
  return str;
}

const std::string& typeName(const uint64_t*) {
  static std::string str("uint64");
  return str;
}

const std::string& typeName(const float*) {
  static std::string str("float");
  return str;
}

const std::string& typeName(const double*) {
  static std::string str("double");
  return str;
}

const std::string& typeName(const std::complex<float>*) {
  static std::string str("fcomplex");
  return str;
}
const std::string& typeName(const std::complex<double>*) {
  static std::string str("dcomplex");
  return str;
}

}  // namespace common
}  // namespace dp3
