// KVpair.h: Implements a KV pair as a pair<string, string>.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_KVPAIR_H
#define LOFAR_COMMON_KVPAIR_H

// Never #include <config.h> or #include <lofar_config.h> in a header file!
#include <ctime>
#include <utility>
#include <iosfwd>
#include <string>

namespace dp3 {
namespace common {

/// \ingroup Common
/// \brief Implements a KV pair as a pair<string, string>.

/// @{

/// Implements a KV pair as a pair<string, string>.
class KVpair : public std::pair<std::string, std::string> {
 public:
  // Note: while this class is not PVSS specific, it is mostly (only?)
  // used by the PVSSGateway, which uses valueType to map the enum values
  // below to PVSS types to query (write).
  // If you add types at all without PVSS support, document that below,
  // so that PVSS users can avoid them.
  KVpair(const std::string& aKey, const std::string& aValue,
         bool genTimestamp = false, bool timestampInKeyname = false);
  KVpair(const std::string& aKey, const char* aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);
  KVpair(const std::string& aKey, bool aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);
  KVpair(const std::string& aKey, int aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);
  KVpair(const std::string& aKey, double aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);
  KVpair(const std::string& aKey, float aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);
  KVpair(const std::string& aKey, time_t aValue, bool genTimestamp = false,
         bool timestampInKeyname = false);

  KVpair();
  ~KVpair();

  // Copying is allowed
  KVpair(const KVpair& that);
  KVpair& operator=(const KVpair& that);
  inline bool operator==(const KVpair& that) const {
    return (first == that.first && second == that.second &&
            timestamp == that.timestamp && valueType == that.valueType);
  }

  // data-members
  double timestamp;  // store also as double
  int16_t valueType;

  enum {
    VT_UNKNOWN = 0,
    VT_STRING,
    VT_BOOL,
    VT_INT,
    VT_DOUBLE,
    VT_FLOAT,
    VT_TIME_T
  };
};

/// @} addgroup

std::ostream& operator<<(std::ostream& os, const KVpair& kv);

}  // namespace common
}  // namespace dp3

#endif
