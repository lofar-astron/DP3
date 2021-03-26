// SubtractMixed.h: Subtract visibilities from buffer after weighting by mixing
// coefficients.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// Subtract visibilities from a buffer after weighting by mixing coefficients.

#ifndef DPPP_SUBTRACTMIXED_H
#define DPPP_SUBTRACTMIXED_H

#include "Baseline.h"
#include "Cursor.h"

namespace dp3 {
namespace base {

/// Subtract visibilities from a buffer after weighting by mixing coefficients.
///
/// \param[in]   nBaseline
/// Number of baselines.
/// \param[in]   nChannel
/// Number of frequency channels.
/// \param[in]   baselines
/// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
/// \param[in]   data
/// A cursor for a 3-D buffer of observed visibilities of shape
/// (\p nBaseline, \p nChannel, 4).
/// \param[in]   model
/// A cursor for a 3-D buffer of simulated visibilities of shape
/// (\p nBaseline, \p nChannel, 4).
/// \param[in]   weight
/// A cursor for a 3-D buffer of mixing weight of shape
/// (\p nBaseline, \p nChannel, 4).
void subtract(size_t nBaseline, size_t nChannel,
              const_cursor<Baseline> baselines, cursor<fcomplex> data,
              const_cursor<dcomplex> model, const_cursor<dcomplex> weight);

}  // namespace base
}  // namespace dp3

#endif
