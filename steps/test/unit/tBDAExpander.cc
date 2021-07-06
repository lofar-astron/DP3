// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../../base/BDABuffer.h"
#include "../../../common/ParameterSet.h"
#include "../../BDAExpander.h"
#include "../../InputStep.h"
#include "mock/MockInput.h"
#include "mock/MockStep.h"

#include <boost/optional.hpp>

#include <string>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::BDAExpander;

const unsigned int kNCorr = 4;
const unsigned int kNChan = 8;
const unsigned int kStartChan = 0;
const unsigned int kNTime = 3;
const double kStartTime = 0.0;
const double kInterval = 2.0;
const double kUVW[3]{1.0, 1.0, 1.0};
const std::size_t kNIntervals = 3;

const int kNAntennas = 3;
const int kNBaselines = 2;
const std::string kMsName{};
const std::string kAntennaSet{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
const std::vector<casacore::MPosition> kAntPos{casacore::MVPosition{0, 0, 0},
                                               casacore::MVPosition{10, 0, 0},
                                               casacore::MVPosition{20, 0, 0}};
const std::vector<double> kAntDiam(kNAntennas, 1.0);
const std::vector<int> kAnt1_2Bl{0, 0};
const std::vector<int> kAnt2_2Bl{1, 2};

BOOST_AUTO_TEST_SUITE(bda_expander, *boost::unit_test::tolerance(0.001) *
                                        boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(time_expansion) {
  std::size_t pool_size = kNCorr * kNChan * kNBaselines * kNIntervals;

  std::unique_ptr<BDABuffer> buffer{new BDABuffer(pool_size)};

  const double bda_first_time = kStartTime + kInterval / 2;
  std::vector<std::size_t> baseline_id{0, 1};

  const std::vector<std::complex<float>> kData1(kNCorr * kNChan,
                                                std::complex<float>(1.0, 1.0));
  const std::vector<std::complex<float>> kData2(kNCorr * kNChan,
                                                std::complex<float>(2.0, 2.0));
  const std::vector<std::complex<float>> kData3(kNCorr * kNChan,
                                                std::complex<float>(3.0, 3.0));
  const std::vector<std::complex<float>> kData4(kNCorr * kNChan,
                                                std::complex<float>(4.0, 4.0));

  const std::vector<float> kWeight(kNCorr * kNChan, 1.0);

  // baseline 0, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[0], kNCorr,
                 kNChan, kData1.data(), nullptr, kWeight.data(), nullptr, kUVW);
  // baseline 0, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan, kData2.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);
  // baseline 0, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan, kData3.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);
  // baseline 1, timeslot = 1 + 2 + 3
  buffer->AddRow(kStartTime + (kNIntervals * kInterval) / 2,
                 kNIntervals * kInterval, kNIntervals * kInterval,
                 baseline_id[1], kNCorr, kNChan, kData4.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);

  DPInfo info;
  std::vector<std::vector<double>> chan_freqs(kNBaselines);
  std::vector<std::vector<double>> chan_widths(kNBaselines);
  for (int k = 0; k < kNBaselines; k++) {
    for (std::size_t i = 0; i < kNChan; i++) {
      chan_freqs[k].push_back(i * 10000.0);
      chan_widths[k].push_back(10000.0);
    }
  }

  info.init(kNCorr, kStartChan, kNChan, kNTime, bda_first_time - kInterval / 2,
            kInterval, kMsName, kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(true);
  info.set(std::move(chan_freqs), std::move(chan_widths));

  dp3::common::ParameterSet parset;
  dp3::steps::MockInput mock_input;
  BDAExpander expander("bdaexpander");

  BOOST_CHECK_NO_THROW(expander.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  expander.setNextStep(mock_step);
  expander.process(std::move(buffer));

  // check the size of the output
  BOOST_REQUIRE_EQUAL(mock_step->GetRegularBuffers().size(), kNIntervals);

  // check if start times are properly set
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0].getTime(),
                    bda_first_time);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1].getTime(),
                    bda_first_time + kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2].getTime(),
                    bda_first_time + 2 * kInterval);

  // check that interval size is same for all times, and equal to minimal size
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0].getExposure(), kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1].getExposure(), kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2].getExposure(), kInterval);

  // check if data is correctly copied
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0].getData()(0, 0, 0),
                    kData1.front());
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1].getData()(0, 0, 0),
                    kData2.front());
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2].getData()(0, 0, 0),
                    kData3.front());
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0].getData()(0, 0, 1),
                    kData4.front());
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1].getData()(0, 0, 1),
                    kData4.front());
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2].getData()(0, 0, 1),
                    kData4.front());
}

BOOST_AUTO_TEST_CASE(frequency_expansion) {
  std::size_t pool_size = kNCorr * kNChan * kNBaselines * kNIntervals;

  std::unique_ptr<BDABuffer> buffer{new BDABuffer(pool_size)};

  const double bda_first_time = kStartTime + kInterval / 2;
  const std::vector<std::size_t> baseline_id{0, 1};

  const std::vector<std::complex<float>> kData1(kNCorr * kNChan / 2,
                                                std::complex<float>(1.0, 1.0));
  const std::vector<std::complex<float>> kData2(kNCorr * kNChan,
                                                std::complex<float>(2.0, 2.0));
  const std::vector<std::complex<float>> kData3(kNCorr * kNChan / 2,
                                                std::complex<float>(3.0, 3.0));
  const std::vector<std::complex<float>> kData4(kNCorr * kNChan,
                                                std::complex<float>(4.0, 4.0));
  const std::vector<std::complex<float>> kData5(kNCorr * kNChan / 2,
                                                std::complex<float>(5.0, 5.0));
  const std::vector<std::complex<float>> kData6(kNCorr * kNChan,
                                                std::complex<float>(6.0, 6.0));

  const std::vector<float> kWeight(kNCorr * kNChan, 1.0);

  // baseline 0, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[0], kNCorr,
                 kNChan / 2, kData1.data(), nullptr, kWeight.data(), nullptr,
                 kUVW);
  // baseline 1, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[1], kNCorr,
                 kNChan, kData2.data(), nullptr, kWeight.data(), nullptr, kUVW);
  // baseline 0, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan / 2, kData3.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);
  // baseline 1, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[1], kNCorr, kNChan, kData4.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);
  // baseline 0, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan / 2, kData5.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);
  // baseline 1, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[1], kNCorr, kNChan, kData6.data(), nullptr,
                 kWeight.data(), nullptr, kUVW);

  DPInfo info;

  std::vector<std::vector<double>> chan_freqs{
      // Baseline 0: channels are averaged
      {10000.0, 30000.0, 50000.0, 70000.0},
      // Baseline 1: channels are not averaged
      {5000.0, 15000.0, 25000.0, 35000.0, 45000.0, 55000.0, 65000.0, 75000.0}};

  std::vector<std::vector<double>> chan_widths{
      // Baseline 0: channels are averaged
      {2 * 10000.0, 2 * 10000.0, 2 * 10000.0, 2 * 10000.0},
      // Baseline 1: channels are not averaged
      {10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0}};

  info.init(kNCorr, kStartChan, kNChan, kNTime, bda_first_time - kInterval / 2,
            kInterval, kMsName, kAntennaSet);

  info.set(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(true);
  info.set(std::move(chan_freqs), std::move(chan_widths));

  dp3::common::ParameterSet parset;
  dp3::steps::MockInput mock_input;
  BDAExpander expander("bdaexpander");

  BOOST_CHECK_NO_THROW(expander.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  expander.setNextStep(mock_step);
  expander.process(std::move(buffer));

  BOOST_REQUIRE_EQUAL(mock_step->GetRegularBuffers().size(), 3u);

  for (auto it = mock_step->GetRegularBuffers().begin();
       it != mock_step->GetRegularBuffers().end(); ++it) {
    BOOST_REQUIRE_EQUAL(it->getData().shape(),
                        casacore::IPosition(3, kNCorr, kNChan, kNBaselines));
  }

  for (unsigned int idx_channel = 0; idx_channel < kNChan; ++idx_channel) {
    for (unsigned int idx_corr = 0; idx_corr < kNCorr; ++idx_corr) {
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[0].getData()(idx_corr, idx_channel, 0),
          kData1.front());
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[0].getData()(idx_corr, idx_channel, 1),
          kData2.front());
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[1].getData()(idx_corr, idx_channel, 0),
          kData3.front());
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[1].getData()(idx_corr, idx_channel, 1),
          kData4.front());
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[2].getData()(idx_corr, idx_channel, 0),
          kData5.front());
      BOOST_CHECK_EQUAL(
          mock_step->GetRegularBuffers()[2].getData()(idx_corr, idx_channel, 1),
          kData6.front());
    }
  }
}

BOOST_AUTO_TEST_CASE(wrong_input_parameters) {
  DPInfo info;
  info.init(kNCorr, kStartChan, kNChan, kNTime, kStartTime, kInterval, kMsName,
            kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(false);

  dp3::common::ParameterSet parset;
  dp3::steps::MockInput mock_input;
  BDAExpander expander("bdaexpander");
  BOOST_CHECK_THROW(expander.updateInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
