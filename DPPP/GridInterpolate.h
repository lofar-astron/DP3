//# GridInterpolate.h: Interpolate data from regular 2d grid to another
//# Copyright (C) 2018
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
//# $Id: GridInterpolate.h 37169 2017-04-19 12:41:21Z dijkema $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_GRIDINTERPOLATE_H
#define DPPP_GRIDINTERPOLATE_H

// @file
// @brief Interpolate data from regular 2d grid to another

#include <vector>
#include <stdexcept>

namespace DP3 {
  //! Get the nearest-neighbor indices
  ///*! \param ax_src[in]    Vector with points where the data is defined.
  //                         Should be increasing.
  // *  \param ax_tgt[in]    Vector with the points at which the values are
  //                         needed.  Should be increasing.
  //  *  \param[out] indices Vector (same length as ax_tgt) with for each number
  //                         in ax_src, the index of the nearest point in ax_src.
  //   *  \param[in] nearest Get the nearest point. If false, gets the largest
  //                         point that is smaller.
  //    */
  void getAxisIndices(const std::vector<double>& ax_src,
                      const std::vector<double>& ax_tgt,
                      std::vector<size_t>& indices,
                      bool nearest = true);

  //! Regrid 2d-gridded data onto another 2d grid
  /*! \param[in] x_src x-axis on which the data is defined
   *  \param[in] y_src y-axis on which the data is defined
   *  \param[in] x_tgt x-axis on which the data will be evaluated
   *  \param[in] y_tgt y-axis on which the data will be evaluated
   *  \param[in] vals_src original data, y-axis varies fastest
   *  \param[out] vals_tgt regridded data, y-axis varies fastest
   */
  void gridNearestNeighbor(const std::vector<double>& x_src,
                           const std::vector<double>& y_src,
                           const std::vector<double>& x_tgt,
                           const std::vector<double>& y_tgt,
                           const double* vals_src,
                           double* vals_tgt,
                           bool nearest = true);
}

#endif
