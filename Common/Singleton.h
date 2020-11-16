// Singleton.h: Implementation of a Meyers singleton class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_SINGLETON_H
#define LOFAR_COMMON_SINGLETON_H

namespace DP3 {
/// \ingroup Common
/// \brief Implementation of a Meyers singleton class.

/// @{

/// Singleton implements the so-called Meyers singleton (see Item 26 in
/// <em>More Effective C++</em>, by Scott Meyers).
///
/// \attention The Meyers singleton is \e not thread-safe. So, you can only
/// safely use this Singleton class <em>as long as</em> there is only one
/// thread running. Note that, in general, the static initialization and
/// destruction phases (i.e. before main() has started and after main() has
/// finished), are single threaded.
///
/// \attention The order of destruction of static objects is
/// undetermined. This Singleton class therefore suffers from the so-called
/// "Dead Reference Problem"
template <typename T>
class Singleton {
 public:
  /// Return a reference to the object \c T. The object \c T is created when
  /// instance() is called for the first time.
  static T& instance();

 private:
  /// @{
  /// Do not allow construction, destruction, copy construction, and
  /// assignment by a third party.
  Singleton();
  Singleton(const Singleton<T>&);
  Singleton<T>& operator=(const Singleton<T>&);
  ~Singleton();
  /// @}
};

/// @}

template <typename T>
T& Singleton<T>::instance() {
  static T obj;
  return obj;
}

}  // namespace DP3

#endif
