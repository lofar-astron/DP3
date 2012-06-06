//# SubtractMixed.h: Subtract visibilities from buffer after weighting by mixing
//# coefficients.
//#
//# Copyright (C) 2012
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#ifndef DPPP_SUBTRACTMIXED_H
#define DPPP_SUBTRACTMIXED_H

// \file
// Subtract visibilities from a buffer after weighting by mixing coefficients.

#include <DPPP/Baseline.h>
#include <DPPP/Cursor.h>
#include <Common/lofar_complex.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{
// Subtract visibilities from a buffer after weighting by mixing coefficients.
//
// \param[in]   nBaseline
// Number of baselines.
// \param[in]   nChannel
// Number of frequency channels.
// \param[in]   baselines
// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
// \param[in]   data
// A cursor for a 3-D buffer of observed visibilities of shape
// (\p nBaseline, \p nChannel, 4).
// \param[in]   model
// A cursor for a 3-D buffer of simulated visibilities of shape
// (\p nBaseline, \p nChannel, 4).
// \param[in]   weight
// A cursor for a 3-D buffer of mixing weight of shape
// (\p nBaseline, \p nChannel, 4).
void subtract(size_t nBaseline, size_t nChannel,
    const_cursor<Baseline> baselines, cursor<fcomplex> data,
    const_cursor<dcomplex> model, const_cursor<dcomplex> weight);
// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
