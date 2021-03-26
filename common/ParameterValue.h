// ParameterValue.h: The value of a parameter
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PARAMETERVALUE_H
#define LOFAR_COMMON_PARAMETERVALUE_H

#include "StringTools.h"

namespace dp3 {
namespace common {

// Forward declaration.
class ParameterRecord;

/// \brief The value of a parameter

/// ParameterValue represent a value of a parameter.
/// It can contain a single value, but also a vector of ParameterValues or
/// a ParameterRecord.
///
/// It contains various functions to obtain the value in the format desired.
class ParameterValue {
 public:
  /// Default constructor uses empty string.
  ParameterValue() {}

  /// Create from the given string.
  /// Optionally left and right whitespace will be removed.
  explicit ParameterValue(const std::string& value, bool trim = true);

  /// Expand the string using StringUtil::expandedArrayString.
  ParameterValue expand() const;

  /// Is the value a vector?
  bool isVector() const { return itsValue[0] == '['; }

  /// Is the value a record?
  bool isRecord() const { return itsValue[0] == '{'; }

  /// Get the value string.
  const std::string& get() const { return itsValue; }

  /// Get the parameter value as a vector of ParameterValues.
  std::vector<ParameterValue> getVector() const;

  /// Get the parameter value as a ParameterRecord.
  ParameterRecord getRecord() const;

  /// Get the parameter value in the given type.
  /// @{
  bool getBool() const { return strToBool(itsValue); }
  int getInt() const { return strToInt(itsValue); }
  unsigned int getUint() const { return strToUint(itsValue); }
  int16_t getInt16() const { return strToInt16(itsValue); }
  uint16_t getUint16() const { return strToUint16(itsValue); }
  int32_t getInt32() const { return strToInt32(itsValue); }
  int32_t getUint32() const { return strToUint32(itsValue); }
  int64_t getInt64() const { return strToInt64(itsValue); }
  uint64_t getUint64() const { return strToUint64(itsValue); }
  float getFloat() const { return strToFloat(itsValue); }
  double getDouble() const { return strToDouble(itsValue); }
  std::string getString() const;
  time_t getTime() const { return StringToTime_t(itsValue); }
  std::vector<bool> getBoolVector() const;
  std::vector<int> getIntVector() const;
  std::vector<unsigned int> getUintVector() const;
  std::vector<int16_t> getInt16Vector() const;
  std::vector<uint16_t> getUint16Vector() const;
  std::vector<int32_t> getInt32Vector() const;
  std::vector<uint32_t> getUint32Vector() const;
  std::vector<int64_t> getInt64Vector() const;
  std::vector<uint64_t> getUint64Vector() const;
  std::vector<float> getFloatVector() const;
  std::vector<double> getDoubleVector() const;
  std::vector<std::string> getStringVector() const;
  std::vector<time_t> getTimeVector() const;
  /// @}

  /// Convert the parameter value to the given type using conversion operators.
  /// @{
  operator bool() const { return getBool(); }
  operator int() const { return getInt(); }
  operator unsigned int() const { return getUint(); }
  operator float() const { return getFloat(); }
  operator double() const { return getDouble(); }
  operator std::string() const { return getString(); }
  operator time_t() const { return getTime(); }
  operator std::vector<bool>() const { return getBoolVector(); }
  operator std::vector<int>() const { return getIntVector(); }
  operator std::vector<unsigned int>() const { return getUintVector(); }
  operator std::vector<float>() const { return getFloatVector(); }
  operator std::vector<double>() const { return getDoubleVector(); }
  operator std::vector<std::string>() const { return getStringVector(); }
  operator std::vector<time_t>() const { return getTimeVector(); }
  /// @{

  /// Convert a string to a time.
  static time_t StringToTime_t(const std::string& aString);

  /// Put or get to/from ostream.
  /// @{
  friend std::ostream& operator<<(std::ostream& os,
                                  const ParameterValue& pval) {
    os << pval.itsValue;
    return os;
  }
  friend std::istream& operator>>(std::istream& os, ParameterValue& pval) {
    os >> pval.itsValue;
    return os;
  }
  /// @}

 private:
  /// Split itsValue into individual values using the comma as separator.
  /// It takes the commas in quoted strings or in compound values into
  /// account.
  std::vector<ParameterValue> splitValue(unsigned int st,
                                         unsigned int last) const;

  /// Return the substring with left and right whitespace removed.
  ParameterValue substr(int st, int end) const {
    return ParameterValue(itsValue.substr(st, end - st));
  }

  std::string itsValue;
};

}  // namespace common
}  // namespace dp3

#endif
