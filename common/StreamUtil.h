// StreamUtil.h: useful stream manipulation methods.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_STREAMUTILX_H
#define LOFAR_COMMON_STREAMUTILX_H

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace dp3 {
namespace common {

/// \brief Useful stream manipulation methods.

/// Handles indentation of text lines. Every time an Indent object is
/// constructed, the static member \c level is incremented by one, and every
/// time an Indent object is destructed \c level is decremented by one. To
/// increment the amount of indentation you simply create an Indent
/// object. When this object goes out of scope, the amount of indentation is
/// automagically decremented.
class Indent {
 public:
  /// Constructor. Increments indentation level.
  Indent() { lvl++; }
  /// Destructor. Decrements indentation level.
  ~Indent() { lvl--; }
  /// Return the amount of indentation.
  static unsigned int level() { return lvl; }
  /// Return the token to be printed per indentation level
  static const std::string& token() { return tok; }

 private:
  /// Indentation level.
  static unsigned int lvl;
  /// Token to be printed per indentation level.
  static const std::string tok;
};

/// Print an indentation that depends on the number of Indent objects
/// currently in existence.
inline std::ostream& indent(std::ostream& os) {
  for (unsigned int i = 0; i < Indent::level(); ++i) {
    os << Indent::token();
  }
  return os;
}

/// Declare the functions.
/// Write a std::pair.
template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p);

/// Write any container to the given std::ostream.
template <typename ITER>
void print(std::ostream& os, ITER begin, ITER end, const char* separator = ",",
           const char* prefix = "[", const char* postfix = "]");

/// Print the contents of a vector enclosed in square brackets, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);

/// Print the contents of a set enclosed in square brackets, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& v);

/// Print the contents of a map enclosed in braces, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::map<T, U>& m);

/// Write a std::pair.
template <typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
  os << '<' << p.first << ',' << p.second << '>';
  return os;
}

/// Write any container to the given ostream.
template <typename ITER>
inline void print(std::ostream& os, ITER begin, ITER end, const char* separator,
                  const char* prefix, const char* postfix) {
  os << prefix;
  if (begin != end) {
    os << *begin;
    ++begin;
  }
  for (; begin != end; ++begin) {
    os << separator << *begin;
  }
  os << postfix;
}

/// Print the contents of a vector enclosed in square brackets, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  print(os, v.begin(), v.end(), ",", "[", "]");
  return os;
}

/// Print the contents of a set enclosed in square brackets, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {
  print(os, s.begin(), s.end(), ",", "[", "]");
  return os;
}

/// Print the contents of a map enclosed in braces, using a comma
/// as separator.
/// \note operator<<() must be defined for type \c T.
template <typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const std::map<T, U>& m) {
  print(os, m.begin(), m.end(), ", ", "{", "}");
  return os;
}

}  // namespace common
}  // namespace dp3

#endif
