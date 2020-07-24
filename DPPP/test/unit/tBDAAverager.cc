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

#include "../../BDAAverager.h"
#include "../../BDABuffer.h"
#include "mock/MockStep.h"

#include <boost/make_unique.hpp>

using DP3::DPPP::BDAAverager;
using DP3::DPPP::BDABuffer;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::DPInfo;

namespace {
const unsigned int kNCorr = 4;
const unsigned int kNChan = 3;
const unsigned int kStartChan = 0;
const unsigned int kNTime = 10;
const double kStartTime = 0.0;
const double kInterval = 1.0;
const std::string kMsName{};
const std::string kAntennaSet{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2", "ant3"};
const std::vector<casacore::MPosition> kAntPos{
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{0, 0, 100},
    casacore::MVPosition{0, 0, 200}, casacore::MVPosition{0, 0, 400}};
const std::vector<double> kAntDiam(4, 1.0);

/**
 * Create a buffer with artifical data values.
 * @param time Start time for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguising buffers.
 */
std::unique_ptr<DPBuffer> createBuffer(const double time,
                                       std::size_t n_baselines,
                                       const float base_value) {
  casacore::Cube<casacore::Complex> data(kNCorr, kNChan, n_baselines);
  float data_value = base_value + 42.0f;
  for (casacore::Complex& element : data) {
    element = {data_value, -data_value};
    data_value += 1.0f;
  }

  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), 1.0f);
  casacore::Cube<bool> full_res_flags(kNChan, 1, n_baselines, false);
  casacore::Matrix<double> uvw(3, n_baselines);
  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    uvw(0, bl) = base_value + bl * 10.0 + 0.0;
    uvw(1, bl) = base_value + bl * 10.0 + 1.0;
    uvw(2, bl) = base_value + bl * 10.0 + 2.0;
  }

  auto buffer = boost::make_unique<DPBuffer>();
  buffer->setTime(time);
  buffer->setExposure(kInterval);
  buffer->setData(data);
  buffer->setWeights(weights);
  buffer->setFlags(flags);
  buffer->setFullResFlags(full_res_flags);
  buffer->setUVW(uvw);

  return buffer;
}

void checkData(const std::complex<float>& expected,
               const std::complex<float>& actual) {
  // BDABuffer::TimeIsEqual should also work with data values.
  BOOST_TEST(expected.real() == actual.real());

  BOOST_CHECK(BDABuffer::TimeIsEqual(expected.imag(), actual.imag()));
}

void checkWeight(const float expected, const float actual) {
  // BDABuffer::TimeIsEqual should also work with weight values.
  BOOST_CHECK(BDABuffer::TimeIsEqual(expected, actual));
}

void checkUvw(const double expected, const double actual) {
  // BDABuffer::TimeIsEqual should also work with uvw values.
  BOOST_CHECK(BDABuffer::TimeIsEqual(expected, actual));
}

void checkBuffer(const DPBuffer& expected, const BDABuffer& actual) {
  const std::size_t n_corr = expected.getData().shape()[0];
  const std::size_t n_chan = expected.getData().shape()[1];
  const std::size_t n_bl = expected.getData().shape()[2];

  // The BDA buffer should have a row for each baseline.
  BOOST_CHECK_EQUAL(n_bl, actual.GetRows().size());

  for (std::size_t bl = 0; bl < n_bl; ++bl) {
    const BDABuffer::Row& row = actual.GetRows()[bl];
    BOOST_CHECK(BDABuffer::TimeIsEqual(expected.getTime(), row.time_));
    BOOST_CHECK(BDABuffer::TimeIsEqual(expected.getExposure(), row.interval_));
    // ??? TODO:compare row_nr ???
    BOOST_REQUIRE_EQUAL(bl, row.baseline_nr_);
    BOOST_REQUIRE_EQUAL(n_chan, row.n_channels_);
    BOOST_REQUIRE_EQUAL(n_corr, row.n_correlations_);
    checkUvw(expected.getUVW()(0, bl), row.uvw_[0]);
    checkUvw(expected.getUVW()(1, bl), row.uvw_[1]);
    checkUvw(expected.getUVW()(2, bl), row.uvw_[2]);

    std::complex<float>* row_data = row.data_;
    bool* row_flag = row.flags_;
    float* row_weight = row.weights_;
    bool* row_full_res_flag = row.full_res_flags_;
    BOOST_REQUIRE(row_data);
    BOOST_REQUIRE(row_flag);
    BOOST_REQUIRE(row_weight);
    BOOST_REQUIRE(row_full_res_flag);
    for (std::size_t chan = 0; chan < n_chan; ++chan) {
      for (std::size_t corr = 0; corr < n_corr; ++corr) {
        checkData(expected.getData()(corr, chan, bl), *row_data);
        BOOST_CHECK_EQUAL(expected.getFlags()(corr, chan, bl), *row_flag);
        checkWeight(expected.getWeights()(corr, chan, bl), *row_weight);
        // !!! TODO: add proper full res flags test.
        BOOST_CHECK_EQUAL(false, *row_full_res_flag);
        ++row_data;
        ++row_flag;
        ++row_weight;
        ++row_full_res_flag;
      }
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bda_averager, *boost::unit_test::tolerance(0.01) *
                                        boost::unit_test::tolerance(0.01f))

BOOST_AUTO_TEST_CASE(single_baseline) {
  const std::size_t kNBaselines = 1;
  const std::vector<int> kAnt1{0};
  const std::vector<int> kAnt2{1};

  BDAAverager averager;
  DPInfo info;
  info.init(kNCorr, kStartChan, kNChan, kNTime, kStartTime, kInterval, kMsName,
            kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, kAnt1, kAnt2);

  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::unique_ptr<DPBuffer> buffer0 =
      createBuffer(kStartTime + 0.0 * kInterval, kNBaselines, 0.0);
  std::unique_ptr<DPBuffer> buffer1 =
      createBuffer(kStartTime + 1.0 * kInterval, kNBaselines, 100.0);
  std::unique_ptr<DPBuffer> buffer2 =
      createBuffer(kStartTime + 2.0 * kInterval, kNBaselines, 200.0);

  auto mock_step = std::make_shared<DP3::DPPP::MockStep>();
  averager.setNextStep(mock_step);

  // When the BDAAverager merely copies data, it may delay output until
  // the next process() call.
  BOOST_CHECK(averager.process(*buffer0));
  BOOST_CHECK(averager.process(*buffer1));
  BOOST_CHECK(mock_step->GetBdaBuffers().size() >= 1);
  BOOST_CHECK(averager.process(*buffer2));
  BOOST_CHECK(mock_step->GetBdaBuffers().size() >= 2);

  mock_step->CheckFinishCount(0);
  averager.finish();
  mock_step->CheckFinishCount(1);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(3));
  checkBuffer(*buffer0, *mock_step->GetBdaBuffers()[0]);
  checkBuffer(*buffer1, *mock_step->GetBdaBuffers()[1]);
  checkBuffer(*buffer2, *mock_step->GetBdaBuffers()[2]);
}

BOOST_AUTO_TEST_SUITE_END()
