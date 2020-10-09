// Copyright (C) 2020
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

#include <boost/test/unit_test.hpp>

#include "../../SolutionInterval.h"
#include "../../DPBuffer.h"
#include "../../MSReader.h"
#include "../../../Common/Timer.h"

#include <boost/make_unique.hpp>
#include <boost/optional.hpp>

using DP3::NSTimer;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::MSReader;
using DP3::DPPP::SolutionInterval;

namespace {
const int kNBL = 2;
const int kNCorr = 1;
const int kNChan = 1;

DPBuffer InitBuffer() {
  DPBuffer buffer;

  // Set uvw
  casacore::Matrix<double> uvw(3, kNBL);
  for (int i = 0; i < kNBL; ++i) {
    uvw(0, i) = i * kNBL + 1;
    uvw(1, i) = i * kNBL + 2;
    uvw(2, i) = i * kNBL + 3;
  }
  buffer.setUVW(uvw);

  // Set data
  casacore::Cube<casacore::Complex> data(kNCorr, kNChan, kNBL);
  for (int i = 0; i < int(data.size()); ++i) {
    data.data()[i] =
        casacore::Complex(i + i * kNBL * 10, i - 1000 + i * kNBL * 6);
  }
  buffer.setData(data);

  // // Set flags
  casacore::Cube<bool> flags(data.shape());
  flags = false;
  buffer.setFlags(flags);

  // // Set FullResFlags
  casacore::Cube<bool> fullResFlags(kNChan, 1, kNBL);
  fullResFlags = false;
  buffer.setFullResFlags(fullResFlags);

  // Set weights
  casacore::Cube<float> weights(kNCorr, kNChan, kNBL);
  weights = 1.;
  buffer.setWeights(weights);

  return buffer;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(solutioninterval)

/// Test if buffer inserted is the same
BOOST_AUTO_TEST_CASE(insertion) {
  MSReader input;
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  size_t n_dirs = 3;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(&input, n_solution, buffer_size, n_dirs, timer);
  solInt.CopyBuffer(buffer);

  BOOST_TEST(&solInt[0] != &buffer);
  BOOST_TEST(solInt.NSolution() == n_solution);
  BOOST_TEST(solInt[0].getData().tovector() == buffer.getData().tovector());
  BOOST_TEST(solInt[0].getFlags().tovector() == buffer.getFlags().tovector());
  BOOST_TEST(solInt[0].getFullResFlags().tovector() ==
             buffer.getFullResFlags().tovector());
  BOOST_TEST(solInt[0].getUVW().tovector() == buffer.getUVW().tovector());
  BOOST_TEST(solInt[0].getWeights().tovector() ==
             buffer.getWeights().tovector());
}

/// Test that the limit cannot be exceeded
BOOST_AUTO_TEST_CASE(limit) {
  MSReader input;
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  size_t n_dirs = 3;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(&input, n_solution, buffer_size, n_dirs, timer);

  solInt.CopyBuffer(buffer);
  BOOST_CHECK_THROW(solInt.CopyBuffer(buffer), std::runtime_error);
}

/// Test if buffer is a copy and can be changed
BOOST_AUTO_TEST_CASE(copy) {
  MSReader input;
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  size_t n_dirs = 3;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(&input, n_solution, buffer_size, n_dirs, timer);
  solInt.CopyBuffer(buffer);

  BOOST_TEST(&solInt[0] != &buffer);
  BOOST_TEST(solInt[0].getData().tovector() == buffer.getData().tovector());
  BOOST_TEST(solInt[0].getFlags().tovector() == buffer.getFlags().tovector());
  BOOST_TEST(solInt[0].getFullResFlags().tovector() ==
             buffer.getFullResFlags().tovector());
  BOOST_TEST(solInt[0].getUVW().tovector() == buffer.getUVW().tovector());
  BOOST_TEST(solInt[0].getWeights().tovector() ==
             buffer.getWeights().tovector());
}

/// Copy a buffer, change a weight and test if it is restored
BOOST_AUTO_TEST_CASE(restore) {
  MSReader input;
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  size_t n_dirs = 3;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(&input, n_solution, buffer_size, n_dirs, timer);
  solInt.CopyBuffer(buffer);

  // Overwrite some values in the buffer
  casacore::Complex new_data(42.0f, -42.0f);
  *solInt.DataPtrs()[0] = new_data;
  float new_weight = 0.5;
  *solInt.WeightPtrs()[0] = new_weight;

  BOOST_TEST(solInt[0].getData().tovector() != buffer.getData().tovector());
  BOOST_TEST(solInt[0].getWeights().tovector() !=
             buffer.getWeights().tovector());

  solInt.RestoreFlagsAndWeights();

  BOOST_TEST(solInt[0].getData().tovector() != buffer.getData().tovector());
  BOOST_TEST(solInt[0].getWeights().tovector() ==
             buffer.getWeights().tovector());
}

/// Test if Fit resized the arrays
BOOST_AUTO_TEST_CASE(fit) {
  MSReader input;
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 2;
  size_t n_dirs = 3;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(&input, n_solution, buffer_size, n_dirs, timer);
  BOOST_TEST(solInt.Size() == 0U);
  solInt.CopyBuffer(buffer);
  solInt.Fit();

  BOOST_TEST(solInt.Size() == 1U);
  BOOST_TEST(solInt.ModelData().size() == 1U);
  BOOST_TEST(solInt.ModelData().capacity() == 2U);
  BOOST_TEST(solInt.DataPtrs().size() == 1U);
  BOOST_TEST(solInt.DataPtrs().capacity() == 2U);
  BOOST_TEST(solInt.ModelDataPtrs().size() == 1U);
  BOOST_TEST(solInt.ModelDataPtrs().capacity() == 2U);
  BOOST_TEST(solInt.WeightPtrs().size() == 1U);
  BOOST_TEST(solInt.WeightPtrs().capacity() == 2U);
}

BOOST_AUTO_TEST_SUITE_END()
