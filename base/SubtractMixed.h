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
/// \param[out] var_before variance before subtraction
/// \param[out] var_after variance after subtraction
void subtract(size_t nBaseline, size_t nChannel,
              const_cursor<Baseline> baselines,
              cursor<std::complex<float>> data,
              const_cursor<std::complex<double>> model,
              const_cursor<std::complex<double>> weight, float &var_before,
              float &var_after);

}  // namespace base
}  // namespace dp3

#endif
