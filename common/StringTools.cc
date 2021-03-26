// StringUtil.cc: implementation of the string utilities class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!
#include "StringTools.h"

#include <boost/algorithm/string.hpp>

#include <cstring>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <cctype>
#include <errno.h>
#include <map>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace dp3 {
namespace common {

enum rangeElementEnum { PARENTHESIS, ASTERISK, DOTDOT };

typedef std::map<std::string, rangeElementEnum> rangeElementTable;
typedef rangeElementTable::const_iterator rangeElementLookup;

std::vector<std::string> stringtools::tokenize(const std::string& str,
                                               const std::string& delims) {
  std::vector<std::string> tokens;
  std::string::size_type pos = 0;
  std::string::size_type pos0;

  while ((pos0 = str.find_first_not_of(delims, pos)) != std::string::npos) {
    pos = str.find_first_of(delims, pos0 + 1);
    if (pos - pos0 > 0) {  // If pos == std::string::npos then substr() clamps.
      tokens.push_back(str.substr(pos0, pos - pos0));
    }
  }

  return tokens;
}

// formatString(format, ...) --> string up to 10Kb
//
// Function that accepts printf like arguments and returns a string.
const std::string formatString(const char* format, ...) {
  char tmp_cstring[10240];
  va_list ap;

  va_start(ap, format);
  vsnprintf(tmp_cstring, sizeof(tmp_cstring), format, ap);
  va_end(ap);

  return std::string(tmp_cstring);
}

unsigned int lskipws(const std::string& value, unsigned int st,
                     unsigned int end) {
  for (; st < end && isspace(value[st]); ++st)
    ;

  return st;
}

unsigned int rskipws(const std::string& value, unsigned int st,
                     unsigned int end) {
  for (; end > st && isspace(value[end - 1]); --end)
    ;

  return end;
}

bool strToBool(const std::string& aString) {
  char firstChar = aString.c_str()[0];
  if ((firstChar == 't') || (firstChar == 'T') || (firstChar == '1') ||
      (firstChar == 'Y') || (firstChar == 'y'))
    return (true);

  if ((firstChar == 'f') || (firstChar == 'F') || (firstChar == '0') ||
      (firstChar == 'N') || (firstChar == 'n'))
    return (false);

  throw std::runtime_error(aString + " is not a boolean value");
}

long strToLong(const std::string& aString) {
  const char* str = aString.c_str();
  int st = lskipws(aString, 0, aString.size());
  int end = rskipws(aString, st, aString.size());
  char* endPtr;
  long val;
  // Clear errno since strtol does not do it.
  errno = 0;
  // We don't want octal values, so use base 10 unless a hex value is given.
  if (end > st + 2 && str[st] == '0' &&
      (str[st + 1] == 'x' || str[st + 1] == 'X')) {
    val = strtol(str + st, &endPtr, 0);
  } else {
    val = strtol(str + st, &endPtr, 10);
  }
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not an integer value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(aString + " is invalid or outside long int range");
  return val;
}

int strToInt(const std::string& aString) {
  long val = strToLong(aString);
  if (sizeof(int) != sizeof(long)) {
    if (sizeof(int) == 4) {
      if (val < long(-32768) * 65536 || val > 2147483647L)
        throw std::runtime_error(std::to_string(val) +
                                 " is outside 4-byte int range");
    }
  }
  return val;
}

int32_t strToInt32(const std::string& aString) {
  long val = strToLong(aString);
  if (sizeof(int32_t) != sizeof(long)) {
    if (val < long(-32768) * 65536 || val > 2147483647L)
      throw std::runtime_error(std::to_string(val) +
                               " is outside int32_t range");
  }
  return val;
}

int16_t strToInt16(const std::string& aString) {
  long val = strToLong(aString);
  if (val < -32768L || val > 32767L)
    throw std::runtime_error(std::to_string(val) + " is outside int16_t range");
  return val;
}

unsigned long strToUlong(const std::string& aString) {
  const char* str = aString.c_str();
  int st = lskipws(aString, 0, aString.size());
  int end = rskipws(aString, st, aString.size());
  char* endPtr;
  unsigned long val;
  // For negative numbers, strtoul does not raise an error:
  // It uses unsigned integer wraparound rules.
  if (str[st] == '-') {
    throw std::runtime_error(aString + " is not a positive value");
  }
  // Clear errno since strtoul does not do it.
  errno = 0;
  // We don't want octal values, so use base 10 unless a hex value is given.
  if (end > st + 2 && str[st] == '0' &&
      (str[st + 1] == 'x' || str[st + 1] == 'X')) {
    val = strtoul(str + st, &endPtr, 0);
  } else {
    val = strtoul(str + st, &endPtr, 10);
  }
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not an integer value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(aString +
                             " is invalid or outside long unsigned int range");
  return val;
}

unsigned int strToUint(const std::string& aString) {
  unsigned long val = strToUlong(aString);
  if (sizeof(unsigned int) != sizeof(unsigned long)) {
    if (sizeof(unsigned int) == 4) {
      if (val > 4294967295UL)
        throw std::runtime_error(std::to_string(val) +
                                 " is outside 4-byte unsigned int range");
    }
  }
  return val;
}

uint32_t strToUint32(const std::string& aString) {
  unsigned long val = strToUlong(aString);
  if (sizeof(uint32_t) != sizeof(unsigned long)) {
    if (val > 4294967295UL)
      throw std::runtime_error(std::to_string(val) +
                               " is outside uint32_t range");
  }
  return val;
}

uint16_t strToUint16(const std::string& aString) {
  unsigned long val = strToUlong(aString);
  if (val > 65535UL)
    throw std::runtime_error(std::to_string(val) +
                             " is outside uint16_t range");
  return val;
}

float strToFloat(const std::string& aString) {
  const char* str = aString.c_str();
  int end = rskipws(aString, 0, aString.size());
  char* endPtr;
  // Clear errno since strtof does not do it.
  errno = 0;
  double val = strtof(str, &endPtr);
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not a floating point value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(aString + " is invalid or outside float range");
  return val;
}

double strToDouble(const std::string& aString) {
  const char* str = aString.c_str();
  int end = rskipws(aString, 0, aString.size());
  char* endPtr;
  // Clear errno since strtod does not do it.
  errno = 0;
  double val = strtod(str, &endPtr);
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not a floating point value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(aString + " is invalid or outside double range");
  return val;
}

static_assert(sizeof(int64_t) == sizeof(long) ||
                  sizeof(int64_t) != sizeof(long long),
              "strToInt64: sizeof(int64) cannot be handled");

int64_t strToInt64(const std::string& aString) {
  if (sizeof(int64_t) == sizeof(long)) return strToLong(aString);
  const char* str = aString.c_str();
  int st = lskipws(aString, 0, aString.size());
  int end = rskipws(aString, st, aString.size());
  char* endPtr;
  long long val;
  // Clear errno since strtoll does not do it.
  errno = 0;
  // We don't want octal values, so use base 10 unless a hex value is given.
  if (end > st + 2 && str[st] == '0' &&
      (str[st + 1] == 'x' || str[st + 1] == 'X')) {
    val = strtoll(str + st, &endPtr, 0);
  } else {
    val = strtoll(str + st, &endPtr, 10);
  }
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not an integer value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(aString +
                             " is invalid or outside long long int range");
  return val;
}

static_assert(sizeof(int64_t) == sizeof(long) ||
                  sizeof(int64_t) != sizeof(long long),
              "strToUint64: sizeof(uint64) cannot be handled");

uint64_t strToUint64(const std::string& aString) {
  if (sizeof(uint64_t) == sizeof(unsigned long)) return strToUlong(aString);
  const char* str = aString.c_str();
  int st = lskipws(aString, 0, aString.size());
  int end = rskipws(aString, st, aString.size());
  char* endPtr;
  unsigned long long val;
  // For negative numbers, strtoul does not raise an error:
  // It uses unsigned integer wraparound rules.
  if (str[st] == '-') {
    throw std::runtime_error(aString + " is not a positive value");
  }
  // Clear errno since strtoull does not do it.
  errno = 0;
  // We don't want octal values, so use base 10 unless a hex value is given.
  if (end > st + 2 && str[st] == '0' &&
      (str[st + 1] == 'x' || str[st + 1] == 'X')) {
    val = strtoull(str + st, &endPtr, 0);
  } else {
    val = strtoull(str + st, &endPtr, 10);
  }
  if (endPtr != str + end)
    throw std::runtime_error(aString + " is not an integer value");
  if (errno == ERANGE || errno == EINVAL)
    throw std::runtime_error(
        aString + " is invalid or outside long long unsigned int range");
  return val;
}

//
// expandedArrayString(string)
//
// Given een array string ( '[ xx..xx, xx ]' ) this utility expands the string
// by replacing ranges with the fill series.
// Eg. [ lii001..lii003, lii005 ] --> [ lii001, lii002, lii003, lii005 ]
//     [ 10*0         ] --> [ 0,0,0,0,0,0,0,0,0,0 ]
//     [ 3*(0;1;2;3)  ] --> [ 0,1,2,3,0,1,2,3,0,1,2,3 ]
//     [ 3*(300..303) ] --> [ 300,301,302,303,300,301,302,303,300,301,302,303 ]
//     [ 2*(5*0)      ] --> [ 0,0,0,0,0,0,0,0,0,0 ]

// ----------------------- ATTENTION !!!----------------------------------
// This routine has been copied to the Navigator software
// (MAC/Navigator/scripts/libs/nav_usr/CS1/CS1_Common.ctl)
// if you change anything structural change the Navigator part also please
// -----------------------------------------------------------------------
std::string expandRangeString(const std::string& strng) {
  std::string str(strng);
  unsigned int i = 0;
  unsigned int last = str.size();
  while (i + 1 < last) {
    if (str[i] == '\'' || str[i] == '"') {
      // Ignore a quoted part.
      std::string::size_type pos = str.find(str[i], i + 1);
      if (pos == std::string::npos)
        throw std::runtime_error("Unbalanced quoted string at position " +
                                 std::to_string(i) + " in " + str);
      i = pos;
    } else if (str[i] == '.' && str[i + 1] == '.') {
      // Found ..; look back for number and prefix.
      // First find number.
      int endnum = rskipws(str, 0, i);
      int stnum = endnum - 1;
      while (stnum >= 0 && isdigit(str[stnum])) {
        --stnum;
      }
      int lennum = endnum - stnum - 1;
      if (lennum > 0) {
        int num = strToInt(str.substr(stnum + 1, lennum));
        // Found number, now find possible prefix.
        // We could say that the prefix has to be alphanumeric, but more
        // general is to accept all characters except ([,*; blank and tab.
        int stalp = stnum;
        while (stalp >= 0 && str[stalp] != '(' && str[stalp] != '[' &&
               str[stalp] != ',' && str[stalp] != ';' && str[stalp] != '*' &&
               str[stalp] != ' ' && str[stalp] != '\t') {
          --stalp;
        }
        stalp++;
        std::string prefix = str.substr(stalp, stnum - stalp + 1);
        std::string suffix;
        // Now find part after the .. which can contain the same prefix.
        // if no parenthesis was used.
        i = lskipws(str, i + 2, last);
        if (prefix.size() > 0 && last - i > prefix.size() &&
            str.substr(i, prefix.size()) == prefix) {
          i += prefix.size();
        }
        // Test if a digit.
        if (isdigit(str[i])) {
          stnum = i;
          // Skip to end of number.
          while (i < last && isdigit(str[i])) {
            ++i;
          }
          endnum = strToInt(str.substr(stnum, i - stnum));
          // We really have something like xxx000..004
          // Find a possible suffix.
          unsigned int stsuf = i;
          while (i < last && str[i] != ')' && str[i] != ']' && str[i] != ',' &&
                 str[i] != ';' && str[i] != '*' && str[i] != ' ' &&
                 str[i] != '\t') {
            ++i;
          }
          if (i > stsuf) {
            suffix = str.substr(stsuf, i - stsuf);
          }
          // Removes braces if the prefix ends and suffix starts with it.
          int lpre = prefix.size();
          int lsuf = suffix.size();
          if (lpre > 0 && lsuf > 0 && prefix[lpre - 1] == '{' &&
              suffix[0] == '}') {
            prefix = prefix.substr(0, lpre - 1);
            suffix = suffix.substr(1, lsuf - 1);
          }
          // Fill it in in ascending or descending order.
          std::ostringstream ostr;
          if (num < endnum) {
            for (; num <= endnum; ++num) {
              ostr << prefix << std::setfill('0') << std::setw(lennum) << num
                   << suffix;
              if (num != endnum) ostr << ',';
            }
          } else {
            for (; num >= endnum; --num) {
              ostr << prefix << std::setfill('0') << std::setw(lennum) << num
                   << suffix;
              if (num != endnum) ostr << ',';
            }
          }
          str.replace(stalp, i - stalp, ostr.str());
          int diff = ostr.str().size() - (i - stalp);
          last += diff;
          // Start again, because something like aa00..01bcd00..03ef might be
          // used.
          i = stalp - 1;
        }
      }
    }
    ++i;
  }
  return str;
}

std::string expandMultString(const std::string& strng) {
  std::string str(strng);
  unsigned int i = 0;
  unsigned int last = str.size();
  while (i + 1 < last) {
    if (str[i] == '\'' || str[i] == '"') {
      // Ignore a quoted part.
      i = skipQuoted(str, i);
    } else {
      if (str[i] != '*') {
        ++i;
      } else {
        // Found *; look back for digits.
        int endnum = rskipws(str, 0, i);
        int stnum = endnum - 1;
        while (stnum >= 0 && isdigit(str[stnum])) {
          --stnum;
        }
        stnum++;
        // Only use it if the number is at begin or preceeded by a delimiter.
        unsigned int j = rskipws(str, 0, stnum);
        int lennum = 0;
        if (j == 0 || str[j - 1] == ',' || str[j - 1] == '[' ||
            str[j - 1] == '(') {
          lennum = endnum - stnum;
        }
        if (lennum == 0) {
          // No number found, so ignore the *.
          ++i;
        } else {
          bool continueAtReplace = true;
          int num = strToInt(str.substr(stnum, lennum));
          // Now find the string that has to be multiplied.
          // This can be a single instance or a set indicated by () or [].
          // First skip possible whitespace after *.
          i = lskipws(str, i + 1, last);
          unsigned int stval = i;
          std::string val;
          if (str[i] == '[') {
            // A set in [].
            i = skipBalanced(str, i, last, ']');
            val = str.substr(stval, i - stval);
          } else if (str[i] == '(') {
            // A set in (). Remove the parentheses.
            i = skipBalanced(str, i, last, ')');
            val = str.substr(stval + 1, i - stval - 2);
            // Replace ; by , (for backward compatibility).
            unsigned int stv = 0;
            while (stv < val.size()) {
              if (val[stv] == '"' || val[stv] == '\'') {
                stv = skipQuoted(val, stv);
              } else {
                if (val[stv] == ';') {
                  val[stv] = ',';
                }
                stv++;
              }
            }
          } else {
            // Any other value is ended by ,]).
            while (i < last) {
              if (str[i] == '"' || str[i] == '\'') {
                // A quoted string.
                i = skipQuoted(str, i);
              } else if (str[i] == '[') {
                // Another string in brackets.
                i = skipBalanced(str, i, last, ']');
              } else if (str[i] == '(') {
                // Another string in parentheses.
                i = skipBalanced(str, i, last, ')');
              } else {
                if (str[i] == ',' || str[i] == ']' || str[i] == ')') {
                  // The end of the value.
                  break;
                }
                // Continue with next character.
                ++i;
              }
            }
            val = str.substr(stval, i - stval);
            ///            continueAtReplace = false;
          }
          // Insert the values num times separated by a comma.
          std::string res;
          res.reserve(num * (val.size() + 1));
          for (int j = 0; j < num; ++j) {
            if (j > 0) res += ',';
            res += val;
          }
          // Replace the value by the new result.
          str.replace(stnum, i - stnum, res);
          last += res.size() - (i - stnum);
          // Continue scanning at start of replace (for possible recursion).
          if (continueAtReplace) {
            i = stnum;
          } else {
            i = stnum + res.size();
          }
        }
      }
    }
  }
  return str;
}

unsigned int skipQuoted(const std::string& str, unsigned int st) {
  std::string::size_type pos = str.find(str[st], st + 1);
  if (pos == std::string::npos)
    throw std::runtime_error("Unbalanced quoted string at position " +
                             std::to_string(st) + " in " + str);
  return pos + 1;
}

unsigned int skipBalanced(const std::string& str, unsigned int st,
                          unsigned int end, char endChar) {
  char ch = str[st++];
  int nrp = 1;
  while (st < end) {
    // First test on end character. In this way it also works well
    // if start and end character are the same.
    if (str[st] == endChar) {
      st++;
      if (--nrp == 0) return st;
    } else if (str[st] == ch) {
      st++;
      nrp++;
    } else if (str[st] == '"' || str[st] == '\'') {
      st = skipQuoted(str, st);
    } else {
      st++;
    }
  }
  throw std::runtime_error(std::string("Unbalanced ") + ch + endChar + " in " +
                           str);
}

std::string expandArrayString(const std::string& str) {
  // Only do it if enclosed in brackets.
  unsigned int st = lskipws(str, 0, str.size());
  unsigned int end = rskipws(str, st, str.size());
  if (st >= end || str[st] != '[' || str[end - 1] != ']') {
    return str;
  }
  // Do expandMult first, otherwise something like 3*lifs001..003 is not
  // handled properly.
  return expandRangeString(expandMultString(str));
}

}  // namespace common
}  // namespace dp3
