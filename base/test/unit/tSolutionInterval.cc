// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../SolutionInterval.h"
#include <dp3/base/DPBuffer.h>
#include "../../../common/Timer.h"
#include "../../../steps/test/unit/mock/MockInput.h"

using dp3::base::DPBuffer;
using dp3::base::SolutionInterval;
using dp3::common::NSTimer;
using dp3::steps::MockInput;

namespace {
const int kNBL = 2;
const int kNCorr = 1;
const int kNChan = 1;

DPBuffer InitBuffer() {
  DPBuffer buffer;

  casacore::Matrix<double> uvw(3, kNBL);
  for (int i = 0; i < kNBL; ++i) {
    uvw(0, i) = i * kNBL + 1;
    uvw(1, i) = i * kNBL + 2;
    uvw(2, i) = i * kNBL + 3;
  }
  buffer.setUVW(uvw);

  casacore::Cube<casacore::Complex> data(kNCorr, kNChan, kNBL);
  for (int i = 0; i < int(data.size()); ++i) {
    data.data()[i] =
        casacore::Complex(i + i * kNBL * 10, i - 1000 + i * kNBL * 6);
  }
  buffer.setData(data);

  casacore::Cube<bool> flags(data.shape());
  flags = false;
  buffer.setFlags(flags);

  casacore::Cube<float> weights(kNCorr, kNChan, kNBL);
  weights = 1.0f;
  buffer.setWeights(weights);

  return buffer;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(solutioninterval)

/// Test if buffer inserted is the same
BOOST_AUTO_TEST_CASE(insertion) {
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(n_solution, buffer_size, timer);
  BOOST_TEST(solInt.Size() == 0U);

  solInt.PushBack(buffer);
  BOOST_TEST_REQUIRE(solInt.Size() == 1U);

  BOOST_TEST(&solInt[0] != &buffer);
  BOOST_TEST(solInt.NSolution() == n_solution);
  BOOST_TEST(solInt[0].GetCasacoreData().tovector() ==
             buffer.GetCasacoreData().tovector());
  BOOST_TEST(solInt[0].GetCasacoreFlags().tovector() ==
             buffer.GetCasacoreFlags().tovector());
  BOOST_TEST(solInt[0].GetCasacoreUvw().tovector() ==
             buffer.GetCasacoreUvw().tovector());
  BOOST_TEST(solInt[0].GetCasacoreWeights().tovector() ==
             buffer.GetCasacoreWeights().tovector());
}

/// Test that the limit cannot be exceeded
BOOST_AUTO_TEST_CASE(limit) {
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(n_solution, buffer_size, timer);

  solInt.PushBack(buffer);
  BOOST_CHECK_THROW(solInt.PushBack(buffer), std::runtime_error);
}

/// Test if buffer is a copy and can be changed
BOOST_AUTO_TEST_CASE(copy) {
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(n_solution, buffer_size, timer);
  solInt.PushBack(buffer);

  BOOST_TEST(&solInt[0] != &buffer);
  BOOST_TEST(solInt[0].GetCasacoreData().tovector() ==
             buffer.GetCasacoreData().tovector());
  BOOST_TEST(solInt[0].GetCasacoreFlags().tovector() ==
             buffer.GetCasacoreFlags().tovector());
  BOOST_TEST(solInt[0].GetCasacoreUvw().tovector() ==
             buffer.GetCasacoreUvw().tovector());
  BOOST_TEST(solInt[0].GetCasacoreWeights().tovector() ==
             buffer.GetCasacoreWeights().tovector());
}

/// Copy a buffer, change a weight and test if it is restored
BOOST_AUTO_TEST_CASE(restore) {
  NSTimer timer;
  size_t n_solution = 0;
  size_t buffer_size = 1;
  DPBuffer buffer = InitBuffer();

  SolutionInterval solInt(n_solution, buffer_size, timer);
  solInt.PushBack(buffer);

  // Overwrite the values in the buffer
  casacore::Complex new_data(42.0f, -42.0f);
  solInt.DataBuffers()[0].GetCasacoreData() = new_data;
  float new_weight = 0.5;
  solInt.DataBuffers()[0].GetCasacoreWeights() = new_weight;

  BOOST_TEST(solInt[0].GetCasacoreData().tovector() !=
             buffer.GetCasacoreData().tovector());
  BOOST_TEST(solInt[0].GetCasacoreWeights().tovector() !=
             buffer.GetCasacoreWeights().tovector());

  solInt.RestoreFlagsAndWeights();

  BOOST_TEST(solInt[0].GetCasacoreData().tovector() !=
             buffer.GetCasacoreData().tovector());
  BOOST_TEST(solInt[0].GetCasacoreWeights().tovector() ==
             buffer.GetCasacoreWeights().tovector());
}

BOOST_AUTO_TEST_SUITE_END()
