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
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
    casacore::MVPosition{200, 0, 0}, casacore::MVPosition{2400, 0, 0}};
const std::vector<double> kAntDiam(4, 1.0);

void InitInfo(BDAAverager& averager, const std::vector<int>& ant1,
              const std::vector<int>& ant2) {
  BOOST_REQUIRE(ant1.size() == ant2.size());

  DPInfo info;
  info.init(kNCorr, kStartChan, kNChan, kNTime, kStartTime, kInterval, kMsName,
            kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, ant1, ant2);

  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
}

/**
 * Create a buffer with artifical data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguising buffers.
 *        For distinguishing baselines, this function adds baseline_nr * 10.0.
 *        When the buffer represents averaged data, the base_value should be
 *        the total of the base values of the original buffers.
 *        This function divides the base_value by the supplied weight so the
 *        caller does not have to do that division.
 * @param weight Weight value for the data values in the buffer.
 */
std::unique_ptr<DPBuffer> CreateBuffer(const double time, const double interval,
                                       std::size_t n_baselines,
                                       const float base_value,
                                       const float weight = 1.0) {
  casacore::Cube<casacore::Complex> data(kNCorr, kNChan, n_baselines);
  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), weight);
  casacore::Cube<bool> full_res_flags(kNChan, 1, n_baselines, false);
  casacore::Matrix<double> uvw(3, n_baselines);
  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    float bl_value = (base_value / weight) + (bl * 10.0);
    for (unsigned int corr = 0; corr < kNCorr; ++corr) {
      for (unsigned int chan = 0; chan < kNChan; ++chan) {
        data(corr, chan, bl) = bl_value;
        bl_value += 1.0f;
      }
    }
    uvw(0, bl) = bl_value + 0.0;
    uvw(1, bl) = bl_value + 1.0;
    uvw(2, bl) = bl_value + 2.0;
  }

  auto buffer = boost::make_unique<DPBuffer>();
  buffer->setTime(time);
  buffer->setExposure(interval);
  buffer->setData(data);
  buffer->setWeights(weights);
  buffer->setFlags(flags);
  buffer->setFullResFlags(full_res_flags);
  buffer->setUVW(uvw);

  return buffer;
}

void CheckData(const std::complex<float>& expected,
               const std::complex<float>& actual) {
  BOOST_TEST(expected.real() == actual.real());
  BOOST_TEST(expected.imag() == actual.imag());
}

/**
 * Validate that a BDA row is correct, by comparing it to an input buffer.
 * @param expected An input buffer. The first baseline is used for the
 *        comparison, regardless of the baseline number inside the row.
 * @param row A row from a BDA output buffer.
 * @param baseline_nr The expected baseline number of the row.
 */
void CheckRow(const DPBuffer& expected, const BDABuffer::Row& row,
              std::size_t baseline_nr) {
  const std::size_t n_corr = expected.getData().shape()[0];
  const std::size_t n_chan = expected.getData().shape()[1];

  BOOST_TEST(expected.getTime() == row.time_);
  BOOST_TEST(expected.getExposure() == row.interval_);
  // ??? TODO:compare row_nr ???
  BOOST_REQUIRE_EQUAL(baseline_nr, row.baseline_nr_);
  BOOST_REQUIRE_EQUAL(n_chan, row.n_channels_);
  BOOST_REQUIRE_EQUAL(n_corr, row.n_correlations_);
  BOOST_TEST(expected.getUVW()(0, 0) == row.uvw_[0]);
  BOOST_TEST(expected.getUVW()(1, 0) == row.uvw_[1]);
  BOOST_TEST(expected.getUVW()(2, 0) == row.uvw_[2]);

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
      CheckData(expected.getData()(corr, chan, 0), *row_data);
      BOOST_TEST(expected.getFlags()(corr, chan, 0) == *row_flag);
      BOOST_TEST(expected.getWeights()(corr, chan, 0) == *row_weight);
      // !!! TODO: add proper full res flags test.
      BOOST_TEST(false == *row_full_res_flag);
      ++row_data;
      ++row_flag;
      ++row_weight;
      ++row_full_res_flag;
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bda_averager, *boost::unit_test::tolerance(0.001) *
                                        boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(single_baseline) {
  // With a single baseline, the averager should not average any values.
  // It only converts regular buffers to BDA buffers and copies all contents.
  const std::size_t kNBaselines = 1;
  const std::vector<int> kAnt1{0};
  const std::vector<int> kAnt2{1};

  BDAAverager averager;
  InitInfo(averager, kAnt1, kAnt2);

  std::unique_ptr<DPBuffer> buffer0 =
      CreateBuffer(kStartTime + 0.0 * kInterval, kInterval, kNBaselines, 0.0);
  std::unique_ptr<DPBuffer> buffer1 =
      CreateBuffer(kStartTime + 1.0 * kInterval, kInterval, kNBaselines, 100.0);
  std::unique_ptr<DPBuffer> buffer2 =
      CreateBuffer(kStartTime + 2.0 * kInterval, kInterval, kNBaselines, 200.0);

  auto mock_step = std::make_shared<DP3::DPPP::MockStep>();
  averager.setNextStep(mock_step);

  // When the BDAAverager merely copies data, each input buffer should
  // generate one output buffer.
  BOOST_TEST(averager.process(*buffer0));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));
  BOOST_TEST(averager.process(*buffer1));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(2));
  BOOST_TEST(averager.process(*buffer2));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(3));

  mock_step->CheckFinishCount(0);
  averager.finish();
  mock_step->CheckFinishCount(1);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(3));
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[0]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[1]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[2]->GetRows().size());
  CheckRow(*buffer0, mock_step->GetBdaBuffers()[0]->GetRows()[0], 0);
  CheckRow(*buffer1, mock_step->GetBdaBuffers()[1]->GetRows()[0], 0);
  CheckRow(*buffer2, mock_step->GetBdaBuffers()[2]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(three_baselines) {
  // This test uses three baselines with lengths of 2400, 300 and 200.
  // The averager should copy the contents of the first, long, baseline.
  // The averager should average the second baseline by a factor of 2.
  // The averager should average the third baseline by a factor of 3.
  const std::size_t kNBaselines = 3;
  const std::vector<int> kAnt1{0, 0, 0};
  const std::vector<int> kAnt2{3, 1, 2};

  BDAAverager averager;
  InitInfo(averager, kAnt1, kAnt2);

  // Create input buffers for the averager. For the first baseline, these
  // buffers also contain the expected output data.
  std::vector<std::unique_ptr<DPBuffer>> buffers;
  for (int i = 0; i < 5; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, i * 100.0));
  }

  // Create the expected output data for the second baseline. The third buffer
  // is written out after the finish() call and is thus not averaged.
  std::vector<std::unique_ptr<DPBuffer>> expected_second;
  expected_second.push_back(CreateBuffer(kStartTime + 0 * kInterval,
                                         2 * kInterval, 1, 10.0 + 110.0, 2.0));
  expected_second.push_back(CreateBuffer(kStartTime + 2 * kInterval,
                                         2 * kInterval, 1, 210.0 + 310.0, 2.0));
  expected_second.push_back(
      CreateBuffer(kStartTime + 4 * kInterval, kInterval, 1, 410.0, 1.0));

  // Create the expected output data for the third baseline.
  // The second buffer is written out after the finish() call and thus
  // contains the average of two input buffers instead of three.
  std::vector<std::unique_ptr<DPBuffer>> expected_third;
  expected_third.push_back(
      CreateBuffer(kStartTime, 3 * kInterval, 1, 20.0 + 120.0 + 220.0, 3.0));
  expected_third.push_back(CreateBuffer(kStartTime + 3 * kInterval,
                                        2 * kInterval, 1, 320.0 + 420.0, 2.0));

  auto mock_step = std::make_shared<DP3::DPPP::MockStep>();
  averager.setNextStep(mock_step);

  // The averager should output BDABuffers with 2 rows, since the BDA buffers
  // then cover about the same interval (on average) as the input buffers.

  // First baseline: 1 row, second: 0 rows, third: 0 rows -> 0.5 BDA buffers.
  BOOST_TEST(averager.process(*buffers[0]));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));

  // First baseline: 2 rows, second: 1 row, third: 0 rows -> 1.5 BDA buffers.
  BOOST_TEST(averager.process(*buffers[1]));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  // First baseline: 3 rows, second: 1 row, third: 1 row -> 2.5 BDA buffers.
  BOOST_TEST(averager.process(*buffers[2]));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(2));

  // First baseline: 4 rows, second: 2 rows, third: 1 row -> 3.5 BDA buffers.
  BOOST_TEST(averager.process(*buffers[3]));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(3));

  // First baseline: 5 rows, second: 2 rows, third: 1 row -> 4 BDA buffers.
  BOOST_TEST(averager.process(*buffers[4]));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(4));

  mock_step->CheckFinishCount(0);
  averager.finish();
  mock_step->CheckFinishCount(1);

  const auto& bdabuffers = mock_step->GetBdaBuffers();
  // First baseline: 5 rows, second: 3 rows, third: 2 rows -> 5 BDA buffers.
  BOOST_REQUIRE_EQUAL(bdabuffers.size(), std::size_t(5));
  for (int i = 0; i < 5; ++i) {
    BOOST_REQUIRE_EQUAL(2u, bdabuffers[i]->GetRows().size());
  }
  CheckRow(*buffers[0], bdabuffers[0]->GetRows()[0], 0);
  CheckRow(*buffers[1], bdabuffers[0]->GetRows()[1], 0);
  CheckRow(*expected_second[0], bdabuffers[1]->GetRows()[0], 1);
  CheckRow(*buffers[2], bdabuffers[1]->GetRows()[1], 0);
  CheckRow(*expected_third[0], bdabuffers[2]->GetRows()[0], 2);
  CheckRow(*buffers[3], bdabuffers[2]->GetRows()[1], 0);
  CheckRow(*expected_second[1], bdabuffers[3]->GetRows()[0], 1);
  CheckRow(*buffers[4], bdabuffers[3]->GetRows()[1], 0);
  CheckRow(*expected_second[2], bdabuffers[4]->GetRows()[0], 1);
  CheckRow(*expected_third[1], bdabuffers[4]->GetRows()[1], 2);
}

BOOST_AUTO_TEST_CASE(baselines_mismatch) {
  // The antenna vectors indicate there is a single baseline.
  const std::vector<int> kAnt1{0};
  const std::vector<int> kAnt2{1};

  BDAAverager averager;
  InitInfo(averager, kAnt1, kAnt2);

  auto mock_step = std::make_shared<DP3::DPPP::MockStep>();
  averager.setNextStep(mock_step);

  // This buffer has three baselines.
  std::unique_ptr<DPBuffer> buffer =
      CreateBuffer(kStartTime, kInterval, 3, 0.0);
  BOOST_CHECK_THROW(averager.process(*buffer), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
