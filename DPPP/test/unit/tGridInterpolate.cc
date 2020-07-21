// tGridInterpolate.cc: test program for GridInterpolate
// Copyright (C) 2010
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
//
// $Id: tGridInterpolate.cc 31423 2015-04-03 14:06:21Z dijkema $
//
// @author Tammo Jan Dijkema

#include <boost/test/unit_test.hpp>

#include "../../GridInterpolate.h"

using std::vector;

BOOST_AUTO_TEST_SUITE(gridinterpolate)

BOOST_AUTO_TEST_CASE( test_gridinterpolate ) {
  vector<double> ax_src = {1,3};
  vector<double> ax_tgt = {0.5, 1.5, 2.5, 3.5};  

  vector<size_t> indices;
  DP3::getAxisIndices(ax_src, ax_tgt, indices);
  BOOST_CHECK_EQUAL(indices.size(), ax_tgt.size());
  BOOST_CHECK_EQUAL(indices[0], size_t {0});
  BOOST_CHECK_EQUAL(indices[1], size_t {0});
  BOOST_CHECK_EQUAL(indices[2], size_t {1});
  BOOST_CHECK_EQUAL(indices[3], size_t {1});

  vector<double> x_src = {2,4,8,10};
  vector<double> y_src = {3,6,12};
  vector<double> x_tgt = {1,3.5,9.5,10};
  vector<double> y_tgt = {4,10};
  vector<double> vals_src(x_src.size()*y_src.size());
  vector<double> vals_tgt(x_tgt.size()*y_tgt.size());
  for (size_t i=0; i<vals_src.size(); ++i) {
    vals_src[i] = i;
  }

  DP3::getAxisIndices(x_src, x_tgt, indices);
  BOOST_CHECK_EQUAL(indices.size(), x_tgt.size());
  BOOST_CHECK_EQUAL(indices[0], size_t {0});
  BOOST_CHECK_EQUAL(indices[1], size_t {1});
  BOOST_CHECK_EQUAL(indices[2], size_t {3});
  BOOST_CHECK_EQUAL(indices[3], size_t {3});

  DP3::gridNearestNeighbor(x_src, y_src, x_tgt, y_tgt, vals_src.data(), vals_tgt.data());

  BOOST_CHECK_EQUAL(vals_tgt[0], vals_src[0]);
  BOOST_CHECK_EQUAL(vals_tgt[1], vals_src[2]);
  BOOST_CHECK_EQUAL(vals_tgt[2], vals_src[3]);
  BOOST_CHECK_EQUAL(vals_tgt[3], vals_src[5]);
  BOOST_CHECK_EQUAL(vals_tgt[4], vals_src[9]);
  BOOST_CHECK_EQUAL(vals_tgt[5],vals_src[11]);
  BOOST_CHECK_EQUAL(vals_tgt[6], vals_src[9]);
  BOOST_CHECK_EQUAL(vals_tgt[7], vals_src[11]);
}

BOOST_AUTO_TEST_SUITE_END()
