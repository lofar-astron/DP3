// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../BDAAverager.h"
#include "../../../base/BDABuffer.h"
#include "../../InputStep.h"
#include "../../../common/ParameterSet.h"
#include "mock/MockStep.h"
#include "mock/MockInput.h"

#include <boost/make_unique.hpp>
#include <boost/optional.hpp>

#include <string>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::BDAAverager;

namespace {
const unsigned int kNCorr = 4;
const unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
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
const std::vector<int> kAnt1_1Bl{0};
const std::vector<int> kAnt2_1Bl{1};
const std::vector<int> kAnt1_3Bl{0, 0, 0};
const std::vector<int> kAnt2_3Bl{3, 1, 2};

void InitParset(dp3::common::ParameterSet& parset,
                boost::optional<double> timebase = boost::none,
                boost::optional<double> frequencybase = boost::none,
                boost::optional<double> maxinterval = boost::none,
                boost::optional<int> minchannels = boost::none) {
  if (timebase) {
    parset.add("timebase", std::to_string(*timebase));
  }
  if (frequencybase) {
    parset.add("frequencybase", std::to_string(*frequencybase));
  }
  if (maxinterval) {
    parset.add("maxinterval", std::to_string(*maxinterval));
  }
  if (minchannels) {
    parset.add("minchannels", std::to_string(*minchannels));
  }
}

void InitInfo(DPInfo& info, const std::vector<int>& ant1,
              const std::vector<int>& ant2, std::size_t n_chan = kNChan) {
  BOOST_REQUIRE_EQUAL(ant1.size(), ant2.size());
  std::vector<double> chan_freqs(n_chan);
  std::vector<double> chan_widths(n_chan, 5000.0);
  for (std::size_t i = 0; i < n_chan; i++) {
    chan_freqs[i] = i * 10000.0;
  }

  info.init(kNCorr, kStartChan, n_chan, kNTime, kStartTime, kInterval, kMsName,
            kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, ant1, ant2);
  info.set(std::move(chan_freqs), std::move(chan_widths));
}

void CheckInfo(
    const DPInfo& info, const std::vector<std::vector<double>>& chan_freqs,
    const std::vector<std::vector<double>>& chan_widths,
    const std::vector<unsigned int>& time_avg = std::vector<unsigned int>()) {
  BOOST_TEST(info.needVisData());
  BOOST_TEST_REQUIRE(info.hasBDAChannels());
  BOOST_TEST_REQUIRE(info.nbaselines() == chan_freqs.size());
  if (!time_avg.empty()) {
    BOOST_TEST_REQUIRE(info.nbaselines() == time_avg.size());
  }
  for (std::size_t bl = 0; bl < info.nbaselines(); ++bl) {
    BOOST_TEST(info.chanFreqs(bl) == chan_freqs[bl]);
    BOOST_TEST(info.chanWidths(bl) == chan_widths[bl]);
    if (!time_avg.empty()) {
      BOOST_TEST(info.ntimeAvg(bl) == time_avg[bl]);
    }
  }
}

void Finish(BDAAverager& averager, dp3::steps::MockStep& mock_step) {
  BOOST_TEST(mock_step.FinishCount() == std::size_t(0));
  averager.finish();
  BOOST_TEST(mock_step.FinishCount() == std::size_t(1));
}

/**
 * Create a buffer with artifical data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguising buffers.
 *        For distinguishing baselines, this function adds baseline_nr * 100.0.
 *        When the buffer represents averaged data, the base_value should be
 *        the total of the base values of the original buffers.
 *        This function divides the base_value by the supplied weight so the
 *        caller does not have to do that division.
 * @param channel_counts List for generating channel data.
 *        For input buffers, this list should contain a 1 for each channel.
 *        When generating expected output data, this list should contain the
 *        number of averaged input buffers for each output buffer.
 * @param weight Weight value for the data values in the buffer.
 */
std::unique_ptr<DPBuffer> CreateBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight = 1.0) {
  casacore::Cube<casacore::Complex> data(kNCorr, channel_counts.size(),
                                         n_baselines);
  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), weight);
  casacore::Cube<bool> full_res_flags(channel_counts.size(), 1, n_baselines,
                                      false);
  casacore::Matrix<double> uvw(3, n_baselines);
  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    // Base value for this baseline.
    const float bl_value = (bl * 100.0) + (base_value / weight);

    std::size_t chan = 0;
    float chan_value = bl_value;  // Base value for a group of channels.
    for (std::size_t ch_count : channel_counts) {
      // For each channel, increase chan_base by 10.0.
      // When ch_count == 1, 'value' should equal chan_base.
      // When ch_count > 1, 'value' should be the average for multiple channels.
      const float value = chan_value + 5.0 * (ch_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        data(corr, chan, bl) = value + corr;
        weights(corr, chan, bl) *= ch_count;
      }
      ++chan;
      chan_value += ch_count * 10.0;
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

std::unique_ptr<DPBuffer> CreateSimpleBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight) {
  casacore::Cube<casacore::Complex> data(kNCorr, channel_counts.size(),
                                         n_baselines);
  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), weight);
  casacore::Cube<bool> full_res_flags(channel_counts.size(), 1, n_baselines,
                                      false);
  casacore::Matrix<double> uvw(3, n_baselines);

  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    for (std::size_t chan = 0; chan < channel_counts[bl]; ++chan) {
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        data(corr, chan, bl) = base_value;
        weights(corr, chan, bl) = weight;
      }
    }
    uvw(0, bl) = 0.0;
    uvw(1, bl) = 1.0;
    uvw(2, bl) = 2.0;
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

  BOOST_TEST(expected.getTime() == row.time);
  BOOST_TEST(expected.getExposure() == row.interval);
  // ??? TODO:compare row_nr ???
  BOOST_REQUIRE_EQUAL(baseline_nr, row.baseline_nr);
  BOOST_REQUIRE_EQUAL(n_chan * n_corr, row.GetDataSize());
  BOOST_TEST(expected.getUVW()(0, 0) == row.uvw[0]);
  BOOST_TEST(expected.getUVW()(1, 0) == row.uvw[1]);
  BOOST_TEST(expected.getUVW()(2, 0) == row.uvw[2]);

  std::complex<float>* row_data = row.data;
  bool* row_flag = row.flags;
  float* row_weight = row.weights;
  bool* row_full_res_flag = row.full_res_flags;
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

dp3::steps::MockInput mock_input;

}  // namespace

BOOST_AUTO_TEST_SUITE(bda_averager, *boost::unit_test::tolerance(0.001) *
                                        boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(set_averaging_params_invalid) {
  dp3::common::ParameterSet parset;
  BDAAverager averager(mock_input, parset, "");

  std::vector<unsigned int> time_avg{2};
  std::vector<std::vector<double>> freqs{{0, 15000, 35000}};
  std::vector<std::vector<double>> widths{{5000, 10000, 10000}};

  BOOST_REQUIRE_THROW(
      averager.set_averaging_params(std::vector<unsigned int>(), freqs, widths),
      std::invalid_argument);
  BOOST_REQUIRE_THROW(averager.set_averaging_params(
                          time_avg, std::vector<std::vector<double>>(), widths),
                      std::invalid_argument);
  BOOST_REQUIRE_THROW(averager.set_averaging_params(
                          time_avg, freqs, std::vector<std::vector<double>>()),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_averaging_params_unsupported) {
  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);

  dp3::common::ParameterSet parset;
  BDAAverager averager(mock_input, parset, "");

  std::vector<unsigned int> time_avg{2};
  // Define an unsupported frequency averaging scheme (2-2-1)
  std::vector<std::vector<double>> freqs{{5000, 25000, 40000}};
  std::vector<std::vector<double>> widths{{10000, 10000, 5000}};

  BOOST_REQUIRE_NO_THROW(
      averager.set_averaging_params(time_avg, freqs, widths));
  BOOST_REQUIRE_THROW(averager.updateInfo(info), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(set_averaging_params) {
  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);

  dp3::common::ParameterSet parset;
  BDAAverager averager(mock_input, parset, "");

  std::vector<unsigned int> time_avg{2};
  // Define a supported frequency averaging scheme (1-2-2)
  std::vector<std::vector<double>> freqs{{0, 15000, 35000}};
  std::vector<std::vector<double>> widths{{5000, 10000, 10000}};

  BOOST_REQUIRE_NO_THROW(
      averager.set_averaging_params(time_avg, freqs, widths));
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), freqs, widths, time_avg);
}

BOOST_AUTO_TEST_CASE(no_averaging) {
  const std::size_t kNBaselines = 1;
  const std::size_t kTimeSteps = 7;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  // With the default options, the averager performs no averaging: It only
  // copies data from DPBuffers into BDABuffers.
  dp3::common::ParameterSet parset;
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {info.chanFreqs()}, {info.chanWidths()});

  std::vector<std::unique_ptr<DPBuffer>> buffers;
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, kChannelCounts, i * 1000.0));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // When the BDAAverager merely copies data, each input buffer should
  // generate one output buffer.
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(averager.process(*buffers[i]));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == (i + 1));
  }

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[i]->GetRows().size());
    CheckRow(*buffers[i], mock_step->GetBdaBuffers()[i]->GetRows()[0], 0);
  }
}

BOOST_AUTO_TEST_CASE(time_averaging) {
  const std::size_t kFactor = 2;  // Averaging factor for this test.
  const std::size_t kNBaselines = 1;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];
  const double averaged_interval = 2 * kInterval;
  const double averaged_start_time = kStartTime - kInterval / 2;
  const double averaged_centroid_time =
      averaged_start_time + averaged_interval / 2;

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {info.chanFreqs()}, {info.chanWidths()});

  std::unique_ptr<DPBuffer> buffer0 =
      CreateBuffer(kStartTime + 0.0 * kInterval, kInterval, kNBaselines,
                   kChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> buffer1 =
      CreateBuffer(kStartTime + 1.0 * kInterval, kInterval, kNBaselines,
                   kChannelCounts, 1000.0);
  std::unique_ptr<DPBuffer> buffer2 =
      CreateBuffer(kStartTime + 2.0 * kInterval, kInterval, kNBaselines,
                   kChannelCounts, 2000.0);
  std::unique_ptr<DPBuffer> average01 =
      CreateBuffer(averaged_centroid_time, averaged_interval, kNBaselines,
                   kChannelCounts, 0.0 + 1000.0, 2.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer0));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(*buffer1));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));
  BOOST_TEST(averager.process(*buffer2));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(2));
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[0]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[1]->GetRows().size());
  CheckRow(*average01, mock_step->GetBdaBuffers()[0]->GetRows()[0], 0);
  CheckRow(*buffer2, mock_step->GetBdaBuffers()[1]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(time_averaging_use_weights) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];
  const double averaged_interval = 2 * kInterval;
  const double averaged_start_time = kStartTime - kInterval / 2;
  const double averaged_centroid_time =
      averaged_start_time + averaged_interval / 2;
  const float weight_1 = 5.0;
  const float weight_2 = 2.0;
  const float value_1 = 100.0;
  const float value_2 = 200.0;
  float weight_averaged = weight_1 + weight_2;
  float value_averaged =
      (value_1 * weight_1 + value_2 * weight_2) / weight_averaged;

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeAveragingFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {info.chanFreqs()}, {info.chanWidths()});

  std::unique_ptr<DPBuffer> buffer0 =
      CreateSimpleBuffer(kStartTime + 0.0 * kInterval, kInterval, kNBaselines,
                         kChannelCounts, value_1, weight_1);
  std::unique_ptr<DPBuffer> buffer1 =
      CreateSimpleBuffer(kStartTime + 1.0 * kInterval, kInterval, kNBaselines,
                         kChannelCounts, value_2, weight_2);
  std::unique_ptr<DPBuffer> average01 =
      CreateSimpleBuffer(averaged_centroid_time, averaged_interval, kNBaselines,
                         kChannelCounts, value_averaged, weight_averaged);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer0));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(*buffer1));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetRows().size(), 1u);
  CheckRow(*average01, mock_step->GetBdaBuffers()[0]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(time_averaging_ignore_weights) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];
  const double averaged_interval = 2 * kInterval;
  const double averaged_start_time = kStartTime - kInterval / 2;
  const double averaged_centroid_time =
      averaged_start_time + averaged_interval / 2;

  const float weight_1 = 5.0;
  const float weight_2 = 2.0;
  const float value_1 = 100.0;
  const float value_2 = 200.0;

  // when ignoring the weights, the averaged weight only takes into account of
  // the averaging factor (and not the input weights)
  float weight_averaged = float(kTimeAveragingFactor);
  float value_averaged = (value_1 + value_2) / weight_averaged;

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeAveragingFactor);
  BDAAverager averager(mock_input, parset, "", false);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {info.chanFreqs()}, {info.chanWidths()});

  std::unique_ptr<DPBuffer> buffer0 =
      CreateSimpleBuffer(kStartTime + 0.0 * kInterval, kInterval, kNBaselines,
                         kChannelCounts, value_1, weight_1);
  std::unique_ptr<DPBuffer> buffer1 =
      CreateSimpleBuffer(kStartTime + 1.0 * kInterval, kInterval, kNBaselines,
                         kChannelCounts, value_2, weight_2);
  std::unique_ptr<DPBuffer> average01 =
      CreateSimpleBuffer(averaged_centroid_time, averaged_interval, kNBaselines,
                         kChannelCounts, value_averaged, weight_averaged);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer0));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(*buffer1));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetRows().size(), 1u);
  CheckRow(*average01, mock_step->GetBdaBuffers()[0]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(channel_averaging) {
  const std::size_t kFactor = 3;  // Averaging factor for this test.
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kInputChannelCounts(7, 1);
  const std::vector<std::size_t> kOutputChannelCounts{2, 2, 3};
  const std::vector<double> kOutputFreqs{5000.0, 25000.0, 50000.0};
  const std::vector<double> kOutputWidths{10000.0, 10000.0, 15000.0};

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl, kInputChannelCounts.size());
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  InitParset(parset, boost::none, baseline_length * kFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {kOutputFreqs}, {kOutputWidths});

  std::unique_ptr<DPBuffer> buffer = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kInputChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), std::size_t(1));
  BOOST_TEST(0u == bda_buffers[0]->GetRemainingCapacity());
  BOOST_REQUIRE_EQUAL(1u, bda_buffers[0]->GetRows().size());
  CheckRow(*averaged, bda_buffers[0]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(mixed_averaging) {
  const std::size_t kChannelFactor = 3;
  const std::size_t kTimeFactor = 2;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kInputChannelCounts(7, 1);
  const std::vector<std::size_t> kOutputChannelCounts{2, 2, 3};
  const std::vector<double> kOutputFreqs{5000.0, 25000.0, 50000.0};
  const std::vector<double> kOutputWidths{10000.0, 10000.0, 15000.0};

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl, kInputChannelCounts.size());
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeFactor,
             baseline_length * kChannelFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfo(), {kOutputFreqs}, {kOutputWidths});

  std::unique_ptr<DPBuffer> buffer0 =
      CreateBuffer(kStartTime + 0.0 * kInterval, kInterval, kNBaselines,
                   kInputChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> buffer1 =
      CreateBuffer(kStartTime + 1.0 * kInterval, kInterval, kNBaselines,
                   kInputChannelCounts, 1000.0);
  std::unique_ptr<DPBuffer> buffer2 =
      CreateBuffer(kStartTime + 2.0 * kInterval, kInterval, kNBaselines,
                   kInputChannelCounts, 2000.0);

  const double averaged_interval = 2 * kInterval;
  const double averaged_start_time = kStartTime - kInterval / 2;
  const double averaged_centroid_time_01 =
      averaged_start_time + averaged_interval / 2;

  std::unique_ptr<DPBuffer> average01 =
      CreateBuffer(averaged_centroid_time_01, averaged_interval, kNBaselines,
                   kOutputChannelCounts, 0.0 + 1000.0, 2.0);
  std::unique_ptr<DPBuffer> average2 =
      CreateBuffer(kStartTime + 2.0 * kInterval, kInterval, kNBaselines,
                   kOutputChannelCounts, 2000.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer0));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(*buffer1));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));
  BOOST_TEST(averager.process(*buffer2));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(2));
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[0]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[1]->GetRows().size());
  CheckRow(*average01, mock_step->GetBdaBuffers()[0]->GetRows()[0], 0);
  CheckRow(*average2, mock_step->GetBdaBuffers()[1]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(three_baselines_time_averaging) {
  // This test uses three baselines with lengths of 2400, 300 and 200.
  // The averager should copy the contents of the first, long, baseline.
  // The averager should average the second baseline by a factor of 2.
  // The averager should average the third baseline by a factor of 3.
  const std::size_t kNBaselines = 3;

  DPInfo info;
  InitInfo(info, kAnt1_3Bl, kAnt2_3Bl);
  const double time_threshold = 600.0;  // Time averaging factors become 1,2,3.
  const double chan_threshold = 100.0;  // No channel averaging.

  dp3::common::ParameterSet parset;
  InitParset(parset, time_threshold, chan_threshold);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  // Create input buffers for the averager. For the first baseline, these
  // buffers also contain the expected output data.
  std::vector<std::unique_ptr<DPBuffer>> buffers;
  for (int i = 0; i < 5; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, kChannelCounts, i * 1000.0));
  }

  // Create the expected output data for the second baseline. The third buffer
  // is written out after the finish() call and is thus not averaged.
  std::vector<std::unique_ptr<DPBuffer>> expected_second;

  expected_second.push_back(
      CreateBuffer(kStartTime - kInterval / 2 + 2 * kInterval / 2,
                   2 * kInterval, 1, kChannelCounts, 100.0 + 1100.0, 2.0));
  expected_second.push_back(CreateBuffer(
      kStartTime - kInterval / 2 + 2 * kInterval + 2 * kInterval / 2,
      2 * kInterval, 1, kChannelCounts, 2100.0 + 3100.0, 2.0));
  expected_second.push_back(CreateBuffer(kStartTime + 4 * kInterval, kInterval,
                                         1, kChannelCounts, 4100.0, 1.0));

  // Create the expected output data for the third baseline.
  // The second buffer is written out after the finish() call and thus
  // contains the average of two input buffers instead of three.
  std::vector<std::unique_ptr<DPBuffer>> expected_third;
  expected_third.push_back(CreateBuffer(
      kStartTime - kInterval / 2 + (3 * kInterval) / 2, 3 * kInterval, 1,
      kChannelCounts, 200.0 + 1200.0 + 2200.0, 3.0));
  expected_third.push_back(CreateBuffer(
      kStartTime - kInterval / 2 + 3 * kInterval + 2 * kInterval / 2,
      2 * kInterval, 1, kChannelCounts, 3200.0 + 4200.0, 2.0));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
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

  Finish(averager, *mock_step);

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

BOOST_AUTO_TEST_CASE(three_baselines_channel_averaging) {
  // This test uses three baselines with lengths of 2400, 300 and 200.
  // The averager should copy the contents of the first, long, baseline.
  // The averager should average the second baseline by a factor of 2.
  // The averager should average the third baseline by a factor of 3.
  const std::size_t kNBaselines = 3;
  const std::vector<std::size_t> kInputChannelCounts(9, 1);
  const std::vector<std::size_t> kOutputChannelCounts2{1, 2, 2, 2, 2};
  const std::vector<std::size_t> kOutputChannelCounts3{3, 3, 3};
  const std::size_t kTimeSteps = 5;

  DPInfo info;
  InitInfo(info, kAnt1_3Bl, kAnt2_3Bl, kInputChannelCounts.size());
  const double time_threshold = 100.0;  // No averaging, baselines are longer.
  const double chan_threshold = 600.0;  // Averaging factors become 1, 2, 3.

  dp3::common::ParameterSet parset;
  InitParset(parset, time_threshold, chan_threshold);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  // Create input buffers and expected output data for the averager.
  // For the first baseline, the input buffer holds the expected output data.
  std::vector<std::unique_ptr<DPBuffer>> buffers;
  std::vector<std::unique_ptr<DPBuffer>> expected2;
  std::vector<std::unique_ptr<DPBuffer>> expected3;
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    const double time = kStartTime + i * kInterval;
    buffers.push_back(CreateBuffer(time, kInterval, kNBaselines,
                                   kInputChannelCounts, i * 1000.0));
    expected2.push_back(CreateBuffer(time, kInterval, 1, kOutputChannelCounts2,
                                     i * 1000.0 + 100.0));
    expected3.push_back(CreateBuffer(time, kInterval, 1, kOutputChannelCounts3,
                                     i * 1000.0 + 200.0));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // Since there is no time averaging, the averager should generate a full
  // output buffer for each input buffer.
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(averager.process(*buffers[i]));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == i + 1);
  }

  Finish(averager, *mock_step);

  const auto& bdabuffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bdabuffers.size(), kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(0u == bdabuffers[i]->GetRemainingCapacity());
    const std::vector<BDABuffer::Row> rows = bdabuffers[i]->GetRows();
    BOOST_REQUIRE_EQUAL(3u, rows.size());
    CheckRow(*buffers[i], rows[0], 0);
    CheckRow(*expected2[i], rows[1], 1);
    CheckRow(*expected3[i], rows[2], 2);
  }
}

BOOST_AUTO_TEST_CASE(shape_mismatch) {
  // The antenna vectors indicate there is a single baseline.
  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  dp3::common::ParameterSet parset;
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // This buffer has three baselines.
  std::unique_ptr<DPBuffer> buffer =
      CreateBuffer(kStartTime, kInterval, 3, kChannelCounts, 0.0);
  BOOST_CHECK_THROW(averager.process(*buffer), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(max_interval) {
  const double kMaxInterval = kInterval * 3.5;
  const std::size_t kFactor = std::floor(kMaxInterval);
  const std::size_t kNBaselines = 1;
  const std::size_t kTimeSteps = 7 * kFactor;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  // The threshold implies an averaging factor of 42, however, kMaxInterval
  // should limit the averaging to a factor of 3.
  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * 42.0, boost::none, kMaxInterval);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(
        averager.process(*CreateBuffer(kStartTime + i * kInterval, kInterval,
                                       kNBaselines, kChannelCounts, 0.0)));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == (i + 1) / kFactor);
  }

  Finish(averager, *mock_step);

  const auto& bdabuffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bdabuffers.size(), kTimeSteps / kFactor);
  for (const auto& bdabuffer : bdabuffers) {
    const std::vector<BDABuffer::Row> rows = bdabuffer->GetRows();
    BOOST_REQUIRE_EQUAL(1u, rows.size());
    BOOST_TEST(kInterval * kFactor == rows[0].interval);
  }
}

BOOST_AUTO_TEST_CASE(min_channels) {
  const std::size_t kMinChannels = 3;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kOutputChannelCounts{1, 2, 2};

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  // The threshold implies that all channels should be averaged into one
  // channel, however, we specify a minimum of three channels per baseline.
  InitParset(parset, boost::none, baseline_length * kNChan, boost::none,
             kMinChannels);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::unique_ptr<DPBuffer> buffer =
      CreateBuffer(kStartTime, kInterval, kNBaselines, kChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer));
  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), std::size_t(1));
  CheckRow(*averaged, bda_buffers[0]->GetRows()[0], 0);
}

BOOST_AUTO_TEST_CASE(zero_values_weight) {
  const std::size_t kMinChannels = 3;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kOutputChannelCounts{1, 2, 2};

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  // The threshold implies that all channels should be averaged into one
  // channel, however, we specify a minimum of three channels per baseline.
  InitParset(parset, boost::none, baseline_length * kNChan, boost::none,
             kMinChannels);
  BDAAverager averager(mock_input, parset, "");
  averager.updateInfo(info);

  std::unique_ptr<DPBuffer> buffer = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kChannelCounts, 0.0, 0.0);
  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(*buffer));
  Finish(averager, *mock_step);

  // A weight of 0 should not result in a nan
  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_TEST(!std::isnan(bda_buffers[0]->GetData()[0].imag()));
  BOOST_TEST(bda_buffers[0]->GetData()[0].imag() == 0);
}

BOOST_AUTO_TEST_CASE(force_buffersize) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;
  const int kNInputBuffers = 20;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeAveragingFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::vector<std::unique_ptr<DPBuffer>> input_buffers;
  for (int i = 0; i < kNInputBuffers; i++) {
    input_buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                         kNBaselines, kChannelCounts, 0.0));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // The output sizes must be multiples of kNCorr * kNChan (4 * 5 = 20).
  std::vector<unsigned int> output_sizes{80, 20, 100};
  for (unsigned int output_size : output_sizes) {
    averager.set_next_desired_buffersize(output_size);
  }

  for (int i = 0; i < kNInputBuffers; i++) {
    BOOST_TEST(averager.process(*input_buffers[i]));
  }

  Finish(averager, *mock_step);

  // With the given input, we expect a total of 200 elements in the output. In
  // this test we want to split the 200 elements in 3 buffers of sizes 80, 20
  // and 100
  BOOST_TEST_REQUIRE(mock_step->GetBdaBuffers().size() == output_sizes.size());
  for (unsigned int i = 0; i < output_sizes.size(); i++) {
    BOOST_CHECK_EQUAL(mock_step->GetBdaBuffers()[i]->GetNumberOfElements(),
                      output_sizes[i]);
  }
}

BOOST_AUTO_TEST_CASE(force_buffersize_smaller_than_output) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;
  const int kNInputBuffers = 20;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeAveragingFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::vector<std::unique_ptr<DPBuffer>> input_buffers;
  for (int i = 0; i < kNInputBuffers; i++) {
    input_buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                         kNBaselines, kChannelCounts, 0.0));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // The output sizes must be multiples of kNCorr * kNChan (4 * 5 = 20).
  std::vector<unsigned int> output_sizes{20, 40};
  for (unsigned int output_size : output_sizes) {
    averager.set_next_desired_buffersize(output_size);
  }

  for (int i = 0; i < kNInputBuffers; i++) {
    BOOST_TEST(averager.process(*input_buffers[i]));
  }

  Finish(averager, *mock_step);

  // With the given input, we expect a total of 200 elements in the output. In
  // this test we force the number of elements in the first two output buffers
  // to be 20 and 40. The remaining output buffers will have the default size.
  BOOST_REQUIRE(mock_step->GetBdaBuffers().size() >= output_sizes.size());
  for (unsigned int i = 0; i < output_sizes.size(); i++) {
    BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[i]->GetNumberOfElements(),
                        output_sizes[i]);
  }
}

BOOST_AUTO_TEST_CASE(force_buffersize_bigger_than_output) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;
  const int kNInputBuffers = 20;

  DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  dp3::common::ParameterSet parset;
  InitParset(parset, baseline_length * kTimeAveragingFactor);
  BDAAverager averager(mock_input, parset, "");
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::vector<std::unique_ptr<DPBuffer>> input_buffers;
  for (int i = 0; i < kNInputBuffers; i++) {
    input_buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                         kNBaselines, kChannelCounts, 0.0));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  std::vector<unsigned int> output_sizes{1000};

  for (unsigned int output_size : output_sizes) {
    averager.set_next_desired_buffersize(output_size);
  }

  for (int i = 0; i < kNInputBuffers; i++) {
    BOOST_TEST(averager.process(*input_buffers[i]));
  }

  Finish(averager, *mock_step);

  // With the given input, we expect a total of 200 elements in the output. In
  // this test we check that output has 200 elements in case we force
  // the output buffersize to a value higher than the maximum (all the elements
  // are forced into one single bdabuffer).
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetNumberOfElements(),
                      200u);
}

BOOST_AUTO_TEST_SUITE_END()
