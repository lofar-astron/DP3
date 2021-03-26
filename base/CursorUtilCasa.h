// CursorUtilCasa.h: Helper functions for creating cursors for CASA arrays.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// \brief Helper functions for creating cursors for CASA arrays.

#ifndef DPPP_CURSORUTILCASA_H
#define DPPP_CURSORUTILCASA_H

#include "Cursor.h"

#include <casacore/casa/Arrays/Array.h>

namespace dp3 {
namespace base {

template <typename T>
cursor<T> casa_cursor(casacore::Array<T> &array) {
  return cursor<T>(array.data(), array.ndim(), array.steps().storage());
}

template <typename T>
cursor<T> casa_cursor(casacore::Array<T> &array,
                      const casacore::IPosition &offset) {
  return cursor<T>(&(array(offset)), array.ndim(), array.steps().storage());
}

template <typename T>
const_cursor<T> casa_const_cursor(const casacore::Array<T> &array) {
  return const_cursor<T>(array.data(), array.ndim(), array.steps().storage());
}

template <typename T>
const_cursor<T> casa_const_cursor(const casacore::Array<T> &array,
                                  const casacore::IPosition &offset) {
  return const_cursor<T>(&(array(offset)), array.ndim(),
                         array.steps().storage());
}

}  // namespace base
}  // namespace dp3

#endif
