//# GridInterpolate.cc: Interpolate data from regular 2d grid to another
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
//# $Id: GridInterpolate.cc 37169 2017-04-19 12:41:21Z dijkema $
//#

#include "GridInterpolate.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <stdexcept>

using namespace std;

namespace DP3 {
  void getAxisIndices(const vector<double>& ax_src,
                      const vector<double>& ax_tgt,
                      vector<size_t>& indices,
                      bool nearest) {
    indices.resize(ax_tgt.size()); 
    if (ax_tgt.empty()) {
      return;
    }
    assert(!ax_src.empty());

    double lowmatch, highmatch;

    vector<double>::const_iterator src_val = ax_src.begin();
    vector<double>::const_iterator tgt_val = ax_tgt.begin();
    vector<size_t>::iterator index_val = indices.begin();

    while (tgt_val != ax_tgt.end()) {
      while (*src_val < *tgt_val && src_val != ax_src.end()) {
        src_val++;
      }
      if (src_val == ax_src.begin()) {
        *index_val = src_val - ax_src.begin();
      } else if (src_val == ax_src.end()) {
        *index_val = src_val - ax_src.begin() - 1;
      } else {
        if (nearest) {
          lowmatch = *(src_val-1);
          highmatch = *src_val;

          if (highmatch - *tgt_val < *tgt_val - lowmatch) {
            *index_val = src_val - ax_src.begin();
          } else {
            *index_val = src_val - ax_src.begin() - 1;
          }
        } else {
          *index_val = src_val - ax_src.begin() - 1;
        }
      }
      tgt_val++; index_val++;
    }
  }

  void gridNearestNeighbor(const vector<double>& x_src,
                           const vector<double>& y_src,
                           const vector<double>& x_tgt,
                           const vector<double>& y_tgt,
                           const double* vals_src,
                           double* vals_tgt,
                           bool nearest) {
    vector<size_t> x_indices;
    vector<size_t> y_indices;
    getAxisIndices(x_src, x_tgt, x_indices, nearest);
    getAxisIndices(y_src, y_tgt, y_indices, nearest);

    size_t nx = x_tgt.size();
    size_t ny = y_tgt.size();
    size_t ny_src = y_src.size();
    // y varies fastest

    if (nearest) {
      for (size_t i=0; i<nx; ++i) {
        for (size_t j=0; j<ny; ++j) {
          vals_tgt[i*ny+j] = vals_src[x_indices[i]*ny_src + y_indices[j]];
        }
      }
    } else {
      for (size_t i=0; i<nx; ++i) {
        for (size_t j=0; j<ny; ++j) {
          size_t y0_idx, y1_idx;
          bool interpolate_y = true;
          if (y_tgt[j]<=y_src.front()) {
            y0_idx = 0;
            y1_idx = 0;
            interpolate_y = false;
          } else if (y_tgt[j]>=y_src.back()) {
            y0_idx = y_src.size()-1;
            y1_idx = y_src.size()-1;
            interpolate_y = false;
          } else {
            y0_idx = y_indices[j];
            y1_idx = y_indices[j]+1;
          }

          double f_y0, f_y1;
          if (x_tgt[i]<=x_src.front()) {
            f_y0 = vals_src[y0_idx];
            f_y1 = vals_src[y1_idx];
          }
          else if (x_tgt[i]>=x_src.back()) {
            f_y0 = vals_src[(x_src.size()-1)*ny_src+y0_idx];
            f_y1 = vals_src[(x_src.size()-1)*ny_src+y1_idx];
          }
          else {
            size_t x0_idx = x_indices[i];
            double x0 = x_src[x0_idx];
            double x1 = x_src[x0_idx+1];
            double x = x_tgt[i];
            f_y0 = vals_src[x0_idx*ny_src+y0_idx] + (x-x0)/(x1-x0) *(vals_src[(x0_idx+1)*ny_src+y0_idx]-vals_src[x0_idx*ny_src+y0_idx]);
            f_y1 = vals_src[x0_idx*ny_src+y1_idx] + (x-x0)/(x1-x0) *(vals_src[(x0_idx+1)*ny_src+y1_idx]-vals_src[x0_idx*ny_src+y1_idx]);
          }
          if (interpolate_y) {
            double y0 = y_src[y0_idx];
            double y1 = y_src[y0_idx+1];
            double y = y_tgt[j];
            vals_tgt[i*ny+j] = f_y0 + (y-y0)/(y1-y0)*(f_y1-f_y0);
          } else {
            vals_tgt[i*ny+j] = f_y0;
          }
        }
      }
    }
  }
}
