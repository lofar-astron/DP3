// Apply.h: Apply station Jones matrices to a set of visibilities.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_APPLY_H
#define DPPP_APPLY_H

#include "Baseline.h"
#include "Cursor.h"

#include <complex>

namespace dp3 {
namespace base {

/// \brief Apply station Jones matrices to a set of visibilities.

/// @{
/// Apply station Jones matrices to a set of visibilities.
///
/// \param[in]   nBaseline
/// Number of baselines.
/// \param[in]   nChannel
/// Number of frequency channels.
/// \param[in]   baselines
/// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
/// \param[in]   coeff
/// A cursor for a 2-D buffer of Jones matrix coefficients of shape
/// (No. of stations, 8). Each station index contained in \p baselines should be
/// a valid index for the first axis of \p coeff.
/// \param[in]   data
/// A cursor for a 3-D buffer of visibilities of shape
/// (\p nBaseline, \p nChannel, 4).
void apply(size_t nBaseline, size_t nChannel, const_cursor<Baseline> baselines,
           const_cursor<double> coeff, cursor<std::complex<double>> data);
/// @}

}  // namespace base
}  // namespace dp3

#endif
