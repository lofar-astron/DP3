// StringTools.h: useful string manipulation methods.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_STRINGUTIL_H
#define LOFAR_COMMON_STRINGUTIL_H

#include <string>
#include <vector>

namespace dp3 {
namespace common {

/// Useful string manipulation methods and classes.
namespace stringtools {

/// Tokenize the string \c str using any character in \c delim as a separation
/// character. The result does not contain empty strings; consecutive delimiter
/// occurrences count as a single delimiter. Any delimiter occurrences at the
/// beginning or end of \c str are ignored.
///
/// For example:
/// \code
///    vector<string> tokens = StringUtil::tokenize( " aa\t bb  ", " \t" )
/// \endcode
//
/// would yield the following vector of strings:
/// \verbatim
///    tokens[0] = "aa"
///    tokens[1] = "bb"
/// \endverbatim
std::vector<std::string> tokenize(const std::string& str,
                                  const std::string& delims);

/// \brief Functor to compare two strings.
/// Strings can be compared case sensitive
/// (\c NORMAL) and case insensitive (\c NOCASE).
/// \attention This class does not handle locales properly. It does string
/// comparison the way \c strcmp and \c strcasecmp (or \c stricmp for that
/// matter) do it.
class Compare {
 public:
  /// String comparison mode.
  enum Mode { NORMAL, NOCASE };

  /// Constructor. Initialize the comparison criterion. Default is "normal"
  /// case sensitive comparison.
  Compare(Mode mode = NORMAL) : itsMode(mode) {}

  /// The comparison operator
  bool operator()(const std::string& s1, const std::string& s2) const {
    if (itsMode == NORMAL)
      return s1 < s2;
    else
      return lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(),
                                     nocaseCompare);
  }

 private:
  /// Helper function to do case insensitive comparison of two chars.
  static bool nocaseCompare(char c1, char c2) {
    return toupper(c1) < toupper(c2);
  }

  /// The current comparison mode.
  Mode itsMode;
};

}  // namespace stringtools

//
// formatString(format, ...) --> string up to 10Kb
//
// The function formatString accepts printf-like arguments and returns a
// formatted string. It can be used e.g. in cout constructions:
// cout << formatString("Opening connection with host %%s", hostName);
// In real life this must be %s ofcourse but doxygen need a double %%.
const std::string formatString(const char* format, ...);

// Skip leading whitespace (blanks and horizontal tabs) starting at st.
// It returns the position of the first non-whitespace character.
// It returns end if all whitespace.
// It can be used in combination with rskipws.
unsigned int lskipws(const std::string& value, unsigned int st,
                     unsigned int end);

// Skip trailing whitespace (blanks and horizontal tabs) starting at end.
// It returns the position after the last non-whitespace character, thus
// value.substr(st, end-st) extracts the significant value.
// It returns st if all whitespace.
unsigned int rskipws(const std::string& value, unsigned int st,
                     unsigned int end);

// Skip past a quoted string.
// The quote character is the first character (at position st).
// Usually the quote character is ' or ", but it could be any other character.
// An exception is thrown if no ending quote character is found.
unsigned int skipQuoted(const std::string& str, unsigned int st);

// Skip past the end of a balanced pair of delimiters where nested pairs
// are also skipped. Delimiters in quoted parts are ignored.
// The starting delimiter is the first character in the string (at position st).
// The ending delimiter is given as an argument.
// The function also works fine if starting and ending delimiter are the same.
/// <br>An exception is thrown if the delimiters are not balanced, thus if no
// end delimiter is found before position end.
/// <br>For example, it will skip something like '[[1,2,3],[4,5,6]]'.
unsigned int skipBalanced(const std::string& str, unsigned int st,
                          unsigned int end, char endChar);

/// @name Strict conversion of string to numeric value
// Convert a string to the value of any of the fundamental arithmetic data
// types. It checks if the entire string is used for the conversion.
// An integer value can also be given in hexadecimal format (e.g. 0x123).
// Leading and trailing whitespace is allowed.
// It checks if an integer value does not exceed the data type range.
/// @{
bool strToBool(const std::string& aString);
long strToLong(const std::string& aString);
int strToInt(const std::string& aString);
int32_t strToInt32(const std::string& aString);
int16_t strToInt16(const std::string& aString);
unsigned long strToUlong(const std::string& aString);
unsigned int strToUint(const std::string& aString);
uint32_t strToUint32(const std::string& aString);
uint16_t strToUint16(const std::string& aString);
int64_t strToInt64(const std::string& aString);
uint64_t strToUint64(const std::string& aString);
float strToFloat(const std::string& aString);
double strToDouble(const std::string& aString);
/// @}

/// @name Manipulate strings containing a array specification
// Array specifications are often entered by the user with ranges
// like 3..32,55..58 For converting such a string to a real vector the spec must
// be expanded so that it contains all elements instead of the ranges. Likewise,
// when you present a array to the user you often want to show a spec with the
// ranges instead of all individual elements. See the ParameterSet document for
// a detailed description of the syntax.

/// @{
// Given a string 'xx..xx, xx' this utility expands the string
// by replacing ranges with the fill series.
// Eg. "lii001..003xx, lii005" --> "lii001xx,lii002xx,lii003xx, lii005"
std::string expandRangeString(const std::string&);

// Given a string like '3*str' this utility expands the string
// by replacing the string 3 times.
// Eg. "3*0" --> "0,0,0"
std::string expandMultString(const std::string&);

// Apply expandMultString and expandRangeString (in that order) for an array
// string which must be enclosed in square brackets.
std::string expandArrayString(const std::string&);
/// @}

}  // namespace common
}  // namespace dp3

#endif
