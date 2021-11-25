// EstimateMixed.h: Estimate Jones matrices for several directions
// simultaneously. A separate data stream is used for each direction. The
// mixing coefficients quantify the influence of each direction on each of the
// other directions (including time and frequency smearing).
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// Estimate Jones matrices for several directions simultaneously. A separate
/// data stream is used for each direction. The mixing coefficients quantify the
/// influence of each direction on each of the other directions (including time
/// and frequency smearing).

#ifndef DPPP_ESTIMATEMIXED_H
#define DPPP_ESTIMATEMIXED_H

#include "Baseline.h"
#include "Cursor.h"

#include <vector>

namespace dp3 {
namespace base {

/// Estimate Jones matrices for several directions simultaneously. A separate
/// data stream is used for each direction. The mixing coefficients quantify the
/// influence of each direction on each of the other directions (including time
/// and frequency smearing).
//
/// \param[in]   nDirection
/// Number of directions to estimate Jones matrices for.
/// \param[in]   nStation
/// Number of stations.
/// \param[in]   nBaseline
/// Number of baselines.
/// \param[in]   nChannel
/// Number of frequency channels.
/// \param[in]   data
/// Vector of length \p nDirection of cursors for 3-D buffers of observed
/// visiblity data of shape (\p nBaseline, \p nChannel, 4).
/// \param[in]   model
/// Vector of length \p nDirection of cursors for 3-D buffers of simulated
/// visiblity data of shape (\p nBaseline, \p nChannel, 4).
/// \param[in]   baselines
/// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
/// \param[in]   flag
/// A cursor for a 3-D buffer of observed visibility flags of shape
/// (\p nBaseline, \p nChannel, 4).
/// \param[in]   weight
/// A cursor for a 3-D buffer of observed visibility weights of shape
/// (\p nBaseline, \p nChannel, 4).
/// \param[in]   mix
/// A cursor for a 5-D buffer of mixing weights of shape
/// (\p nBaseline, \p nChannel, 4, \p nDirection, \p nDirection).
/// \param[in]   unknowns
/// A pointer to a buffer of unknowns of size nDirection * nStation * 8.
bool estimate(size_t nDirection, size_t nStation, size_t nBaseline,
              size_t nChannel, const_cursor<Baseline> baselines,
              std::vector<const_cursor<fcomplex>> data,
              std::vector<const_cursor<dcomplex>> model,
              const_cursor<bool> flag, const_cursor<float> weight,
              const_cursor<dcomplex> mix, double* unknowns,
              size_t maxiter = 50);

#ifdef HAVE_LIBDIRAC
/// method for LBFGS (robust) solver, with extra input variables:
///  \param[in] lbfgs_mem
///  LBFGS history size, as a multiple of the size of parameters (unknowns)
///  \param[in] robust_nu
///  Robust noise model degrees of freedom, when >30, it is almost Gaussian
/// The remaining input variables are similar to estimate() method above.
bool estimate(std::size_t n_direction, std::size_t n_station,
              std::size_t n_baseline, std::size_t n_channel,
              const_cursor<Baseline> baselines,
              std::vector<const_cursor<fcomplex>> data,
              std::vector<const_cursor<dcomplex>> model,
              const_cursor<bool> flag, const_cursor<float> weight,
              const_cursor<dcomplex> mix, double* unknowns,
              std::size_t lbfgs_mem, double robust_nu,
              std::size_t max_iter = 50);
#endif /* HAVE_LIBDIRAC */

/// Compute a map that contains the index of the unknowns related to the
/// specified baseline in the list of all unknowns.
/// \param[in] n_direction
/// \param[in] n_station
/// \param[in] baseline
///    A tuple of baseline station ids
/// \param[out] index
///    A vector of size 4*8*n_direction to store the indices
void makeIndex(std::size_t n_direction, std::size_t n_station,
               const Baseline& baseline, unsigned int* index);

}  // namespace base
}  // namespace dp3

#endif
