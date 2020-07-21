// CursorUtilCasa.h: Helper functions for creating cursors for CASA arrays.
//
// Copyright (C) 2012
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// \file
/// \brief Helper functions for creating cursors for CASA arrays.

#ifndef DPPP_CURSORUTILCASA_H
#define DPPP_CURSORUTILCASA_H

#include "Cursor.h"

#include <casacore/casa/Arrays/Array.h>

namespace DP3
{
namespace DPPP
{

template <typename T>
cursor<T> casa_cursor(casacore::Array<T> &array)
{
    return cursor<T>(array.data(), array.ndim(), array.steps().storage());
}

template <typename T>
cursor<T> casa_cursor(casacore::Array<T> &array, const casacore::IPosition &offset)
{
    return cursor<T>(&(array(offset)), array.ndim(), array.steps().storage());
}

template <typename T>
const_cursor<T> casa_const_cursor(const casacore::Array<T> &array)
{
    return const_cursor<T>(array.data(), array.ndim(), array.steps().storage());
}

template <typename T>
const_cursor<T> casa_const_cursor(const casacore::Array<T> &array,
    const casacore::IPosition &offset)
{
    return const_cursor<T>(&(array(offset)), array.ndim(),
        array.steps().storage());
}

} // namespace DPPP
} // namespace LOFAR

#endif
