// GridInterpolate.h: Interpolate data from regular 2d grid to another
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Interpolate data from regular 2d grid to another
/// @author Tammo Jan Dijkema

#ifndef DPPP_GRIDINTERPOLATE_H
#define DPPP_GRIDINTERPOLATE_H

#include <vector>
#include <stdexcept>

namespace DP3 {
/**
 * Get the nearest-neighbor indices
 * \param ax_src[in] Vector with points where the data is defined.
 *                   Should be increasing.
 * \param ax_tgt[in] Vector with the points at which the values are
 *                   needed.  Should be increasing.
 * \param[out] indices Vector (same length as ax_tgt) with for each number
 *                     in ax_src, the index of the nearest point in ax_src.
 * \param[in] nearest Get the nearest point. If false, gets the largest
 *                    point that is smaller.
 */
void getAxisIndices(const std::vector<double>& ax_src,
                    const std::vector<double>& ax_tgt,
                    std::vector<size_t>& indices, bool nearest = true);

/**
 * Regrid 2d-gridded data onto another 2d grid
 * \param[in] x_src x-axis on which the data is defined
 * \param[in] y_src y-axis on which the data is defined
 * \param[in] x_tgt x-axis on which the data will be evaluated
 * \param[in] y_tgt y-axis on which the data will be evaluated
 * \param[in] vals_src original data, y-axis varies fastest
 * \param[out] vals_tgt regridded data, y-axis varies fastest
 */
void gridNearestNeighbor(const std::vector<double>& x_src,
                         const std::vector<double>& y_src,
                         const std::vector<double>& x_tgt,
                         const std::vector<double>& y_tgt,
                         const double* vals_src, double* vals_tgt,
                         bool nearest = true);
}  // namespace DP3

#endif
