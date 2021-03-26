// ParameterValue.h: The value of a parameter
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParameterValue.h"
#include "ParameterRecord.h"

#include <cstdio>
#include <stdexcept>

typedef std::runtime_error ParException;

namespace dp3 {
namespace common {

ParameterValue::ParameterValue(const std::string& value, bool trim)
    : itsValue(value) {
  if (trim) {
    unsigned int st = lskipws(itsValue, 0, itsValue.size());
    unsigned int end = rskipws(itsValue, st, itsValue.size());
    if (st > 0 || end < itsValue.size()) {
      itsValue = itsValue.substr(st, end - st);
    }
  }
}

ParameterValue ParameterValue::expand() const {
  return ParameterValue(expandArrayString(itsValue));
}

std::vector<ParameterValue> ParameterValue::splitValue(
    unsigned int st, unsigned int last) const {
  // Allocate result.
  // Empty result if only whitespace left.
  std::vector<ParameterValue> result;
  st = lskipws(itsValue, st, last);
  if (st == last) {
    return result;
  }
  // Split on commas, but take quotes, braces, parentheses, and brackets
  // into account.
  int nrpar = 0;
  int nrbracket = 0;
  int nrbrace = 0;
  unsigned int i = st;
  while (i < last) {
    if (itsValue[i] == '\'' || itsValue[i] == '"') {
      i = skipQuoted(itsValue, i);
    } else {
      if (itsValue[i] == '(') {
        nrpar++;
      } else if (itsValue[i] == ')') {
        nrpar--;
      } else if (itsValue[i] == '[') {
        if (nrpar != 0)
          throw ParException("Unbalanced () around '" + itsValue + '\'');
        nrbracket++;
      } else if (itsValue[i] == ']') {
        nrbracket--;
      } else if (itsValue[i] == '{') {
        if (nrpar != 0)
          throw ParException("Unbalanced () around '" + itsValue + '\'');
        nrbrace++;
      } else if (itsValue[i] == '}') {
        nrbrace--;
      } else if (itsValue[i] == ',') {
        if (nrpar + nrbracket + nrbrace == 0) {
          result.push_back(ParameterValue(substr(st, i)));
          st = i + 1;
        }
      }
      if (nrpar < 0 || nrbracket < 0 || nrbrace < 0)
        throw ParException("Unbalanced () [] or {} in '" + itsValue + '\'');
      i++;
    }
  }
  result.push_back(ParameterValue(substr(st, last)));
  if (nrpar != 0 || nrbracket != 0 || nrbrace != 0)
    throw ParException("Unbalanced () [] or {} in '" + itsValue + '\'');
  return result;
}

std::vector<ParameterValue> ParameterValue::getVector() const {
  unsigned int st = 1;
  unsigned int last = itsValue.size() - 1;
  // An empty string is an empty vector.
  if (itsValue.empty()) {
    return std::vector<ParameterValue>();
  }
  // A single value if there is no opening and closing bracket.
  if (!(itsValue[0] == '[' && itsValue[last] == ']')) {
    return std::vector<ParameterValue>(1, *this);
  }
  if (itsValue.size() < 2 || itsValue[0] != '[' || itsValue[last] != ']')
    throw ParException("Invalid vector specification in value '" + itsValue +
                       '\'');
  return splitValue(st, last);
}

ParameterRecord ParameterValue::getRecord() const {
  unsigned int st = 1;
  unsigned int last = itsValue.size() - 1;
  if (itsValue.size() < 2 || itsValue[0] != '{' || itsValue[last] != '}')
    throw ParException("Invalid record specification in value '" + itsValue +
                       '\'');
  std::vector<ParameterValue> values(splitValue(st, last));
  // Loop over all values and split in names and values.
  ParameterRecord result;
  for (std::vector<ParameterValue>::const_iterator iter = values.begin();
       iter != values.end(); ++iter) {
    const std::string& str = iter->get();
    unsigned int st = 0;
    if (str[0] == '"' || str[0] == '\'') {
      st = skipQuoted(str, 0);
    }
    std::string::size_type pos = str.find(':', st);
    if (pos == std::string::npos)
      throw ParException("Invalid record specification in value '" + str +
                         '\'');
    // Get name and value.
    // ParameterValue is used to remove possible whitespace.
    // getString also removes quotes in case the name was quoted.
    std::string name(ParameterValue(str.substr(0, pos)).getString());
    std::string value(ParameterValue(str.substr(pos + 1)).get());
    result.add(name, value);
  }
  return result;
}

std::string ParameterValue::getString() const {
  // Remove possible quotes used to escape special chars in the value.
  std::string result;
  unsigned int end = itsValue.size();
  unsigned int pos = 0;
  unsigned int stv = 0;
  while (pos < end) {
    if (itsValue[pos] == '"' || itsValue[pos] == '\'') {
      if (stv < pos) {
        // Add unquoted part.
        result += itsValue.substr(stv, pos - stv);
      }
      // Add quoted part without the quotes.
      stv = pos + 1;
      pos = skipQuoted(itsValue, pos);
      result += itsValue.substr(stv, pos - stv - 1);
      stv = pos;
    } else {
      pos++;
    }
  }
  if (stv < end) {
    // Add remaining part.
    result += itsValue.substr(stv, end - stv);
  }
  return result;
}

std::vector<bool> ParameterValue::getBoolVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<bool> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<int> ParameterValue::getIntVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<int> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<unsigned int> ParameterValue::getUintVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<unsigned int> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<int16_t> ParameterValue::getInt16Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<int16_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getInt16());
  }
  return result;
}

std::vector<uint16_t> ParameterValue::getUint16Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<uint16_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getUint16());
  }
  return result;
}

std::vector<int32_t> ParameterValue::getInt32Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<int32_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getInt32());
  }
  return result;
}

std::vector<uint32_t> ParameterValue::getUint32Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<uint32_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getUint32());
  }
  return result;
}

std::vector<int64_t> ParameterValue::getInt64Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<int64_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getInt64());
  }
  return result;
}

std::vector<uint64_t> ParameterValue::getUint64Vector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<uint64_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(iter->getUint64());
  }
  return result;
}

std::vector<float> ParameterValue::getFloatVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<float> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<double> ParameterValue::getDoubleVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<double> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<std::string> ParameterValue::getStringVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<std::string> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

std::vector<time_t> ParameterValue::getTimeVector() const {
  std::vector<ParameterValue> vec(getVector());
  std::vector<time_t> result;
  result.reserve(vec.size());
  for (std::vector<ParameterValue>::const_iterator iter = vec.begin();
       iter != vec.end(); ++iter) {
    result.push_back(*iter);
  }
  return result;
}

time_t ParameterValue::StringToTime_t(const std::string& aString) {
  time_t theTime;
  char unit[1024];
  unit[0] = '\0';
  if (sscanf(aString.c_str(), "%ld%s", &theTime, unit) < 1) {
    throw ParException(aString + " is not a time value");
  }
  switch (unit[0]) {
    case 's':
    case 'S':
    case '\0':
      break;
    case 'm':
    case 'M':
      theTime *= 60;
      break;
    case 'h':
    case 'H':
      theTime *= 3600;
      break;
  }
  return theTime;
}

}  // namespace common
}  // namespace dp3
