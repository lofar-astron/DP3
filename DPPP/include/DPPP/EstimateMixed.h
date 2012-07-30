//# EstimateMixed.h: Estimate Jones matrices for several directions
//# simultaneously. A separate data stream is used for each direction. The
//# mixing coefficients quantify the influence of each direction on each of the
//# other directions (including time and frequency smearing).
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

#ifndef DPPP_ESTIMATEMIXED_H
#define DPPP_ESTIMATEMIXED_H

// \file
// Estimate Jones matrices for several directions simultaneously. A separate
// data stream is used for each direction. The mixing coefficients quantify the
// influence of each direction on each of the other directions (including time
// and frequency smearing).

#include <DPPP/Baseline.h>
#include <DPPP/Cursor.h>
#include <Common/lofar_complex.h>
#include <Common/lofar_vector.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

// Estimate Jones matrices for several directions simultaneously. A separate
// data stream is used for each direction. The mixing coefficients quantify the
// influence of each direction on each of the other directions (including time
// and frequency smearing).
//
// \param[in]   nDirection
// Number of directions to estimate Jones matrices for.
// \param[in]   nStation
// Number of stations.
// \param[in]   nBaseline
// Number of baselines.
// \param[in]   nChannel
// Number of frequency channels.
// \param[in]   data
// Vector of length \p nDirection of cursors for 3-D buffers of observed
// visiblity data of shape (\p nBaseline, \p nChannel, 4).
// \param[in]   model
// Vector of length \p nDirection of cursors for 3-D buffers of simulated
// visiblity data of shape (\p nBaseline, \p nChannel, 4).
// \param[in]   baselines
// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
// \param[in]   flag
// A cursor for a 3-D buffer of observed visibility flags of shape
// (\p nBaseline, \p nChannel, 4).
// \param[in]   weight
// A cursor for a 3-D buffer of observed visibility weights of shape
// (\p nBaseline, \p nChannel, 4).
// \param[in]   mix
// A cursor for a 5-D buffer of mixing weights of shape
// (\p nBaseline, \p nChannel, 4, \p nDirection, \p nDirection).
// \param[in]   unknowns
// A pointer to a bufer of unknowns of size nDirection * nStation * 8.
// \param[in]   errors
// A pointer to a buffer of errors of size nDirection * nStation * 8. If this
// pointer is null (0), errors will not be estimated.
bool estimate(size_t nDirection, size_t nStation, size_t nBaseline,
    size_t nChannel, const_cursor<Baseline> baselines,
    vector<const_cursor<fcomplex> > data, vector<const_cursor<dcomplex> > model,
    const_cursor<bool> flag, const_cursor<float> weight,
    const_cursor<dcomplex> mix, double *unknowns);
// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
