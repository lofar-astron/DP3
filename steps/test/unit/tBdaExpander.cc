// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <dp3/base/BdaBuffer.h>
#include "../../../common/ParameterSet.h"
#include "../../BdaExpander.h"
#include "mock/MockInput.h"
#include "mock/MockStep.h"

#include <boost/optional.hpp>

#include <string>

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;
using dp3::steps::BdaExpander;

const unsigned int kNCorr = 4;
const unsigned int kNChan = 8;
const unsigned int kNTime = 3;
const double kStartTime = 0.0;
const double kInterval = 2.0;
const double kUVW[3]{1.0, 1.0, 1.0};
const std::size_t kNIntervals = 3;

const int kNAntennas = 3;
const int kNBaselines = 2;
const std::string kMsName{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
const std::vector<casacore::MPosition> kAntPos{casacore::MVPosition{0, 0, 0},
                                               casacore::MVPosition{10, 0, 0},
                                               casacore::MVPosition{20, 0, 0}};
const std::vector<double> kAntDiam(kNAntennas, 1.0);
const std::vector<int> kAnt1_2Bl{0, 0};
const std::vector<int> kAnt2_2Bl{1, 2};

void CheckData(const DPBuffer& buffer,
               const std::vector<std::complex<float>>& data, int baseline) {
  BOOST_CHECK_EQUAL(buffer.GetData().size() % data.size(), 0u);
  for (unsigned int chan = 0; chan < kNChan; ++chan) {
    for (unsigned int corr = 0; corr < kNCorr; ++corr) {
      // When frequency averaging is present in the input data, the size of the
      // data in the bda input buffer is smaller than the size of the expanded
      // regular buffer. For this reason the indexing of the 'data' variable
      // should restart from the beginning once it reaches its size.
      BOOST_CHECK_CLOSE(buffer.GetData()(baseline, chan, corr),
                        data[(chan * kNCorr + corr) % data.size()], 1e-3);
    }
  }
}

void CheckWeights(const std::vector<std::unique_ptr<DPBuffer>>& buffers,
                  const std::vector<std::vector<float>>& weights) {
  float total_weight_after_expansion = 0.0;
  float total_weight_before_expansion = 0.0;
  for (unsigned int i = 0; i < buffers.size(); i++) {
    for (unsigned int baseline = 0; baseline < kNBaselines; ++baseline) {
      for (unsigned int chan = 0; chan < kNChan; ++chan) {
        for (unsigned int corr = 0; corr < kNCorr; ++corr) {
          total_weight_after_expansion +=
              buffers[i]->GetWeights()(baseline, chan, corr);
        }
      }
    }
  }

  for (const std::vector<float>& row_weights : weights) {
    total_weight_before_expansion +=
        std::accumulate(row_weights.begin(), row_weights.end(), 0.0);
  }

  BOOST_CHECK_CLOSE(total_weight_after_expansion, total_weight_before_expansion,
                    1e-3);
}

BOOST_AUTO_TEST_SUITE(bda_expander, *boost::unit_test::tolerance(0.001) *
                                        boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(time_expansion) {
  const std::size_t kPoolSize = kNCorr * kNChan * kNBaselines * kNIntervals;

  auto buffer = std::make_unique<BdaBuffer>(
      kPoolSize,
      Fields(Fields::Single::kData) | Fields(Fields::Single::kWeights));

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

  const std::vector<std::vector<float>> kWeights{
      {std::vector<float>(kNCorr * kNChan, 0.7),
       std::vector<float>(kNCorr * kNChan, 1.5),
       std::vector<float>(kNCorr * kNChan, 3.0),
       std::vector<float>(kNCorr * kNChan, 0.01)}};

  // baseline 0, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[0], kNCorr,
                 kNChan, kData1.data(), nullptr, kWeights[0].data(), kUVW);
  // baseline 0, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan, kData2.data(), nullptr,
                 kWeights[1].data(), kUVW);
  // baseline 0, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan, kData3.data(), nullptr,
                 kWeights[2].data(), kUVW);
  // baseline 1, timeslot = 1 + 2 + 3
  buffer->AddRow(kStartTime + (kNIntervals * kInterval) / 2,
                 kNIntervals * kInterval, kNIntervals * kInterval,
                 baseline_id[1], kNCorr, kNChan, kData4.data(), nullptr,
                 kWeights[3].data(), kUVW);

  std::vector<std::vector<double>> chan_freqs(kNBaselines);
  std::vector<std::vector<double>> chan_widths(kNBaselines);
  for (int k = 0; k < kNBaselines; k++) {
    for (std::size_t i = 0; i < kNChan; i++) {
      chan_freqs[k].push_back(i * 10000.0);
      chan_widths[k].push_back(10000.0);
    }
  }

  DPInfo info(kNCorr, kNChan);
  info.setTimes(bda_first_time, bda_first_time + (kNTime - 1) * kInterval,
                kInterval);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(true);
  info.setChannels(std::move(chan_freqs), std::move(chan_widths));

  dp3::common::ParameterSet parset;
  dp3::steps::MockInput mock_input;
  BdaExpander expander("bdaexpander");

  BOOST_CHECK_NO_THROW(expander.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  expander.setNextStep(mock_step);
  expander.process(std::move(buffer));

  // check the size of the output
  BOOST_REQUIRE_EQUAL(mock_step->GetRegularBuffers().size(), kNIntervals);

  // check if start times are properly set
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0]->GetTime(),
                    bda_first_time);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1]->GetTime(),
                    bda_first_time + kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2]->GetTime(),
                    bda_first_time + 2 * kInterval);

  // check that interval size is same for all times, and equal to minimal size
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0]->GetExposure(),
                    kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[1]->GetExposure(),
                    kInterval);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[2]->GetExposure(),
                    kInterval);

  // check if data is correctly copied
  CheckData(*mock_step->GetRegularBuffers()[0], kData1, 0);
  CheckData(*mock_step->GetRegularBuffers()[1], kData2, 0);
  CheckData(*mock_step->GetRegularBuffers()[2], kData3, 0);
  CheckData(*mock_step->GetRegularBuffers()[0], kData4, 1);
  CheckData(*mock_step->GetRegularBuffers()[1], kData4, 1);
  CheckData(*mock_step->GetRegularBuffers()[2], kData4, 1);

  // check if weights are correctly processed
  CheckWeights(mock_step->GetRegularBuffers(), kWeights);
}

BOOST_AUTO_TEST_CASE(frequency_expansion) {
  std::size_t kPoolSize = kNCorr * kNChan * kNBaselines * kNIntervals;

  auto buffer = std::make_unique<BdaBuffer>(
      kPoolSize,
      Fields(Fields::Single::kData) | Fields(Fields::Single::kWeights));

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

  const std::vector<std::vector<float>> kWeights{
      {std::vector<float>(kNCorr * kNChan / 2, 0.9),
       std::vector<float>(kNCorr * kNChan, 1.3),
       std::vector<float>(kNCorr * kNChan / 2, 9.7),
       std::vector<float>(kNCorr * kNChan, 0.01),
       std::vector<float>(kNCorr * kNChan / 2, 2.0),
       std::vector<float>(kNCorr * kNChan, 4.3)}};

  // baseline 0, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[0], kNCorr,
                 kNChan / 2, kData1.data(), nullptr, kWeights[0].data(), kUVW);
  // baseline 1, timeslot = 1
  buffer->AddRow(bda_first_time, kInterval, kInterval, baseline_id[1], kNCorr,
                 kNChan, kData2.data(), nullptr, kWeights[1].data(), kUVW);
  // baseline 0, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan / 2, kData3.data(), nullptr,
                 kWeights[2].data(), kUVW);
  // baseline 1, timeslot = 2
  buffer->AddRow(bda_first_time + kInterval, kInterval, kInterval,
                 baseline_id[1], kNCorr, kNChan, kData4.data(), nullptr,
                 kWeights[3].data(), kUVW);
  // baseline 0, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[0], kNCorr, kNChan / 2, kData5.data(), nullptr,
                 kWeights[4].data(), kUVW);
  // baseline 1, timeslot = 3
  buffer->AddRow(bda_first_time + 2 * kInterval, kInterval, kInterval,
                 baseline_id[1], kNCorr, kNChan, kData6.data(), nullptr,
                 kWeights[5].data(), kUVW);

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

  DPInfo info(kNCorr, kNChan);
  info.setTimes(bda_first_time, bda_first_time + (kNTime - 1) * kInterval,
                kInterval);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(true);
  info.setChannels(std::move(chan_freqs), std::move(chan_widths));

  dp3::common::ParameterSet parset;
  dp3::steps::MockInput mock_input;
  BdaExpander expander("bdaexpander");

  BOOST_CHECK_NO_THROW(expander.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  expander.setNextStep(mock_step);
  expander.process(std::move(buffer));

  BOOST_REQUIRE_EQUAL(mock_step->GetRegularBuffers().size(), 3u);

  for (const std::unique_ptr<DPBuffer>& buffer :
       mock_step->GetRegularBuffers()) {
    BOOST_REQUIRE_EQUAL(buffer->GetData().shape(0), kNBaselines);
    BOOST_REQUIRE_EQUAL(buffer->GetData().shape(1), kNChan);
    BOOST_REQUIRE_EQUAL(buffer->GetData().shape(2), kNCorr);
  }

  // check if data is correctly copied
  CheckData(*mock_step->GetRegularBuffers()[0], kData1, 0);
  CheckData(*mock_step->GetRegularBuffers()[1], kData3, 0);
  CheckData(*mock_step->GetRegularBuffers()[2], kData5, 0);
  CheckData(*mock_step->GetRegularBuffers()[0], kData2, 1);
  CheckData(*mock_step->GetRegularBuffers()[1], kData4, 1);
  CheckData(*mock_step->GetRegularBuffers()[2], kData6, 1);

  // check if weights are correctly processed
  CheckWeights(mock_step->GetRegularBuffers(), kWeights);
}

BOOST_AUTO_TEST_CASE(wrong_input_parameters) {
  DPInfo info(kNCorr, kNChan);
  info.setTimes(kStartTime, kStartTime + kInterval, kInterval);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1_2Bl, kAnt2_2Bl);
  info.setIsBDAIntervalFactorInteger(false);

  BdaExpander expander("bdaexpander");
  BOOST_CHECK_THROW(expander.updateInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
