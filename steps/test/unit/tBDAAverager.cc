// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../BDAAverager.h"

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <dp3/base/BdaBuffer.h>
#include "../../../common/ParameterSet.h"
#include "mock/MockStep.h"
#include "tStepCommon.h"

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::BdaAverager;
using dp3::steps::test::kTrueFalseRange;

namespace {
const unsigned int kNCorr = 4;
const unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
const unsigned int kStartChan = 0;
const double kStartTime = 0.0;
const double kInterval = 1.0;
const double kFirstTime = 0.5;
const double kLastTime = 9.5;
const std::string kExtraData = "extraData";
const std::string kMsName{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2", "ant3"};
const std::vector<casacore::MPosition> kAntPos{
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
    casacore::MVPosition{200, 0, 0}, casacore::MVPosition{2400, 0, 0}};
const std::vector<double> kAntDiam(4, 1.0);
const std::vector<int> kAnt1_1Bl{0};
const std::vector<int> kAnt2_1Bl{1};
const std::vector<int> kAnt1_3Bl{0, 0, 0};
const std::vector<int> kAnt2_3Bl{3, 1, 2};
const std::string kPrefix = "bda_averager.";

dp3::common::ParameterSet GetParset(std::optional<double> timebase = {},
                                    std::optional<double> frequencybase = {},
                                    std::optional<double> maxinterval = {},
                                    std::optional<int> minchannels = {}) {
  dp3::common::ParameterSet parset;
  parset.add("steps", "[bda_averager]");
  parset.add("bda_averager.type", "bdaaverager");
  if (timebase) {
    parset.add(kPrefix + "timebase", std::to_string(*timebase));
  }
  if (frequencybase) {
    parset.add(kPrefix + "frequencybase", std::to_string(*frequencybase));
  }
  if (maxinterval) {
    parset.add(kPrefix + "maxinterval", std::to_string(*maxinterval));
  }
  if (minchannels) {
    parset.add(kPrefix + "minchannels", std::to_string(*minchannels));
  }
  return parset;
}

DPInfo InitInfo(const std::vector<int>& ant1, const std::vector<int>& ant2,
                std::size_t n_chan = kNChan) {
  BOOST_REQUIRE_EQUAL(ant1.size(), ant2.size());
  std::vector<double> chan_freqs(n_chan);
  std::vector<double> chan_widths(n_chan, 5000.0);
  for (std::size_t i = 0; i < n_chan; i++) {
    chan_freqs[i] = i * 10000.0;
  }

  DPInfo info(kNCorr, kNChan);
  info.setTimes(kFirstTime, kLastTime, kInterval);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, ant1, ant2);
  info.setChannels(std::move(chan_freqs), std::move(chan_widths));
  return info;
}

void CheckInfo(
    const DPInfo& info, const std::vector<std::vector<double>>& chan_freqs,
    const std::vector<std::vector<double>>& chan_widths,
    const std::vector<unsigned int>& time_avg = std::vector<unsigned int>()) {
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

void Finish(BdaAverager& averager, dp3::steps::MockStep& mock_step) {
  BOOST_TEST(mock_step.FinishCount() == std::size_t(0));
  averager.finish();
  BOOST_TEST(mock_step.FinishCount() == std::size_t(1));
}

/**
 * Create a buffer with artificial data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguishing buffers.
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
 *
 * @note The function has been copied from @ref steps/test/unit/tIDGPredict.cc.
 */
std::unique_ptr<dp3::base::DPBuffer> CreateBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight = 1.0) {
  const std::array<std::size_t, 3> kShape{n_baselines, channel_counts.size(),
                                          kNCorr};

  auto buffer = std::make_unique<dp3::base::DPBuffer>(time, interval);
  buffer->GetData().resize(kShape);
  buffer->AddData(kExtraData);
  buffer->GetWeights().resize(kShape);
  buffer->GetFlags().resize(kShape);
  buffer->GetUvw().resize({n_baselines, 3});

  buffer->GetFlags().fill(false);
  buffer->GetWeights().fill(weight);

  for (std::size_t baseline = 0; baseline < n_baselines; ++baseline) {
    // Base value for this baseline.
    const float baseline_value =
        10'000.0 + (baseline * 100.0) + (base_value / weight);

    std::size_t channel = 0;
    float channel_value = baseline_value;  // Base value for a group of channels
    for (std::size_t channel_count : channel_counts) {
      // For each channel, increase channel_value by 10.0.
      // When channel_count == 1, 'value' should equal channel_value.
      // When channel_count > 1, 'value' should be the average for multiple
      // channels.
      const float value = channel_value + 5.0 * (channel_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        buffer->GetData()(baseline, channel, corr) =
            std::complex<float>(value + corr, value + corr + 100.0f);
        buffer->GetData(kExtraData)(baseline, channel, corr) =
            std::complex<float>(2 * (value + corr), value + corr + 200.0f);
        buffer->GetWeights()(baseline, channel, corr) *= channel_count;
      }
      ++channel;
      channel_value += channel_count * 10.0;
    }
    buffer->GetUvw()(baseline, 0) = baseline_value + 0.0;
    buffer->GetUvw()(baseline, 1) = baseline_value + 1.0;
    buffer->GetUvw()(baseline, 2) = baseline_value + 2.0;
  }

  return buffer;
}

std::unique_ptr<DPBuffer> CreateSimpleBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight) {
  auto buffer = std::make_unique<DPBuffer>(time, interval);

  const std::array<std::size_t, 3> kShape{n_baselines, channel_counts.size(),
                                          kNCorr};
  buffer->GetData().resize(kShape);
  buffer->GetWeights().resize(kShape);
  buffer->GetFlags().resize(kShape);
  buffer->GetUvw().resize({n_baselines, 3});

  buffer->GetData().fill(std::complex{0, 0});
  buffer->GetWeights().fill(weight);
  buffer->GetFlags().fill(false);

  DPBuffer::DataType& data = buffer->GetData();
  DPBuffer::WeightsType& weights = buffer->GetWeights();
  DPBuffer::UvwType& uvw = buffer->GetUvw();

  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    for (std::size_t chan = 0; chan < channel_counts[bl]; ++chan) {
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        data(bl, chan, corr) = base_value;
        weights(bl, chan, corr) = weight;
      }
    }
    uvw(bl, 0) = 0.0;
    uvw(bl, 1) = 1.0;
    uvw(bl, 2) = 2.0;
  }

  return buffer;
}

/**
 * Validate that a BDA row is correct, by comparing it to an input buffer.
 * @param expected An input buffer. The first baseline is used for the
 *        comparison, regardless of the baseline number inside the row.
 * @param row A row from a BDA output buffer.
 * @param baseline_nr The expected baseline number of the row.
 */
void CheckRow(const DPBuffer& expected, const BdaBuffer& buffer,
              std::size_t row_index, std::size_t baseline_nr) {
  const BdaBuffer::Row& row = buffer.GetRows()[row_index];

  const std::size_t n_corr = expected.GetWeights().shape(2);
  const std::size_t n_chan = expected.GetWeights().shape(1);

  const std::vector<std::string> data_names = expected.GetDataNames();
  BOOST_TEST(data_names == buffer.GetDataNames());
  BOOST_TEST(expected.GetTime() == row.time);
  BOOST_TEST(expected.GetExposure() == row.interval);
  // ??? TODO:compare row_nr ???
  BOOST_REQUIRE_EQUAL(baseline_nr, row.baseline_nr);
  BOOST_REQUIRE_EQUAL(n_chan * n_corr, row.GetDataSize());
  BOOST_TEST(expected.GetUvw()(0, 0) == row.uvw[0]);
  BOOST_TEST(expected.GetUvw()(0, 1) == row.uvw[1]);
  BOOST_TEST(expected.GetUvw()(0, 2) == row.uvw[2]);

  std::vector<const std::complex<float>*> row_data;
  for (const std::string& name : data_names) {
    row_data.push_back(buffer.GetData(row_index, name));
    BOOST_REQUIRE(row_data.back());
  }
  const bool* row_flag = buffer.GetFlags(row_index);
  const float* row_weight = buffer.GetWeights(row_index);
  BOOST_REQUIRE(row_flag);
  BOOST_REQUIRE(row_weight);

  for (std::size_t chan = 0; chan < n_chan; ++chan) {
    for (std::size_t corr = 0; corr < n_corr; ++corr) {
      for (std::size_t i = 0; i < data_names.size(); ++i) {
        BOOST_TEST(expected.GetData(data_names[i])(0, chan, corr).real() ==
                   row_data[i]->real());
        BOOST_TEST(expected.GetData(data_names[i])(0, chan, corr).imag() ==
                   row_data[i]->imag());
        ++row_data[i];
      }
      BOOST_TEST(expected.GetFlags()(0, chan, corr) == *row_flag);
      BOOST_TEST(expected.GetWeights()(0, chan, corr) == *row_weight);
      ++row_flag;
      ++row_weight;
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bda_averager, *boost::unit_test::tolerance(0.001) *
                                        boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(finish_without_process) {
  const dp3::common::ParameterSet parset = GetParset();
  BdaAverager averager(parset, kPrefix);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  averager.finish();
  BOOST_TEST(mock_step->GetBdaBuffers().empty());
}

BOOST_AUTO_TEST_CASE(required_fields) {
  using dp3::steps::Step;
  dp3::common::ParameterSet parset;

  // By default, BdaAverager uses weights and flags.
  const dp3::common::Fields kExpectedFieldsWeightsFlags =
      Step::kDataField | Step::kFlagsField | Step::kWeightsField |
      Step::kUvwField;
  const BdaAverager averager(parset, kPrefix);
  BOOST_TEST(averager.getRequiredFields() == kExpectedFieldsWeightsFlags);
  // Test with an explicit 'true' argument for the constructor.
  const BdaAverager explicit_weights_flags(parset, kPrefix, true);
  BOOST_TEST(explicit_weights_flags.getRequiredFields() ==
             kExpectedFieldsWeightsFlags);

  // Disabling using weights and flags affects the requirements.
  const BdaAverager no_weights_flags(parset, kPrefix, false);
  BOOST_TEST(no_weights_flags.getRequiredFields() ==
             (Step::kDataField | Step::kUvwField));
}

BOOST_AUTO_TEST_CASE(set_averaging_params_invalid) {
  const dp3::common::ParameterSet parset = GetParset();
  BdaAverager averager(parset, kPrefix);

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
  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);

  const dp3::common::ParameterSet parset = GetParset(2.0);
  BdaAverager averager(parset, "bda_averager.");

  std::vector<unsigned int> time_avg{2};
  // Define an unsupported frequency averaging scheme (2-2-1)
  std::vector<std::vector<double>> freqs{{5000, 25000, 40000}};
  std::vector<std::vector<double>> widths{{10000, 10000, 5000}};

  BOOST_REQUIRE_NO_THROW(
      averager.set_averaging_params(time_avg, freqs, widths));
  BOOST_REQUIRE_THROW(averager.updateInfo(info), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(set_averaging_params) {
  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);

  const dp3::common::ParameterSet parset = GetParset(2.0);
  BdaAverager averager(parset, kPrefix);

  std::vector<unsigned int> time_avg{2};
  // Define a supported frequency averaging scheme (1-2-2)
  std::vector<std::vector<double>> freqs{{0, 15000, 35000}};
  std::vector<std::vector<double>> widths{{5000, 10000, 10000}};

  BOOST_REQUIRE_NO_THROW(
      averager.set_averaging_params(time_avg, freqs, widths));
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), freqs, widths, time_avg);
}

BOOST_AUTO_TEST_CASE(no_averaging) {
  const std::size_t kNBaselines = 1;
  const std::size_t kTimeSteps = 7;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  // With the default options, the averager performs no averaging: It only
  // copies data from DPBuffers into BdaBuffers.
  const dp3::common::ParameterSet parset = GetParset();
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {info.chanFreqs()}, {info.chanWidths()});

  std::vector<std::unique_ptr<DPBuffer>> buffers;
  std::vector<std::unique_ptr<DPBuffer>> expected;
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, kChannelCounts, i * 1000.0));
    expected.push_back(std::make_unique<DPBuffer>(*buffers.back()));
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // When the BdaAverager merely copies data, each input buffer should
  // generate one output buffer.
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(averager.process(std::move(buffers[i])));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == (i + 1));
  }

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[i]->GetRows().size());
    CheckRow(*expected[i], *mock_step->GetBdaBuffers()[i], 0, 0);
  }
}

BOOST_DATA_TEST_CASE(time_averaging, kTrueFalseRange, use_data) {
  const std::size_t kFactor = 2;  // Averaging factor for this test.
  const std::size_t kNBaselines = 1;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];
  const double averaged_interval = 2 * kInterval;
  const double averaged_start_time = kStartTime - kInterval / 2;
  const double averaged_centroid_time =
      averaged_start_time + averaged_interval / 2;

  const dp3::common::ParameterSet parset = GetParset(baseline_length * kFactor);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {info.chanFreqs()}, {info.chanWidths()});

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
  // The buffer average02 is equal to buffer2
  std::unique_ptr<DPBuffer> average02 = std::make_unique<DPBuffer>(*buffer2);

  if (!use_data) {  // Run the test without any visibilities.
    buffer0->RemoveData();
    buffer1->RemoveData();
    buffer2->RemoveData();
    average01->RemoveData();
    average02->RemoveData();
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(std::move(buffer0)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(std::move(buffer1)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));
  BOOST_TEST(averager.process(std::move(buffer2)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(2));
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[0]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[1]->GetRows().size());
  CheckRow(*average01, *mock_step->GetBdaBuffers()[0], 0, 0);
  CheckRow(*average02, *mock_step->GetBdaBuffers()[1], 0, 0);
}

BOOST_AUTO_TEST_CASE(time_averaging_use_weights) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
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

  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * kTimeAveragingFactor);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {info.chanFreqs()}, {info.chanWidths()});

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

  BOOST_TEST(averager.process(std::move(buffer0)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(std::move(buffer1)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetRows().size(), 1u);
  CheckRow(*average01, *mock_step->GetBdaBuffers()[0], 0, 0);
}

BOOST_AUTO_TEST_CASE(time_averaging_ignore_weights) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
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

  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * kTimeAveragingFactor);
  BdaAverager averager(parset, kPrefix, false);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {info.chanFreqs()}, {info.chanWidths()});

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

  BOOST_TEST(averager.process(std::move(buffer0)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(std::move(buffer1)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetRows().size(), 1u);
  CheckRow(*average01, *mock_step->GetBdaBuffers()[0], 0, 0);
}

BOOST_DATA_TEST_CASE(channel_averaging, kTrueFalseRange, use_data) {
  const std::size_t kFactor = 3;  // Averaging factor for this test.
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kInputChannelCounts(7, 1);
  const std::vector<std::size_t> kOutputChannelCounts{2, 2, 3};
  const std::vector<double> kOutputFreqs{5000.0, 25000.0, 50000.0};
  const std::vector<double> kOutputWidths{10000.0, 10000.0, 15000.0};

  const DPInfo info =
      InitInfo(kAnt1_1Bl, kAnt2_1Bl, kInputChannelCounts.size());
  const double baseline_length = info.getBaselineLengths()[0];

  const dp3::common::ParameterSet parset =
      GetParset(std::nullopt, baseline_length * kFactor);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {kOutputFreqs}, {kOutputWidths});

  std::unique_ptr<DPBuffer> buffer = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kInputChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);
  if (!use_data) {  // Run the test without any visibilities.
    buffer->RemoveData();
    averaged->RemoveData();
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(std::move(buffer)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), std::size_t(1));
  BOOST_TEST(0u == bda_buffers[0]->GetRemainingCapacity());
  BOOST_REQUIRE_EQUAL(1u, bda_buffers[0]->GetRows().size());
  CheckRow(*averaged, *bda_buffers[0], 0, 0);
}

BOOST_AUTO_TEST_CASE(mixed_averaging) {
  const std::size_t kChannelFactor = 3;
  const std::size_t kTimeFactor = 2;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kInputChannelCounts(7, 1);
  const std::vector<std::size_t> kOutputChannelCounts{2, 2, 3};
  const std::vector<double> kOutputFreqs{5000.0, 25000.0, 50000.0};
  const std::vector<double> kOutputWidths{10000.0, 10000.0, 15000.0};

  const DPInfo info =
      InitInfo(kAnt1_1Bl, kAnt2_1Bl, kInputChannelCounts.size());
  const double baseline_length = info.getBaselineLengths()[0];

  const dp3::common::ParameterSet parset = GetParset(
      baseline_length * kTimeFactor, baseline_length * kChannelFactor);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));
  CheckInfo(averager.getInfoOut(), {kOutputFreqs}, {kOutputWidths});

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

  BOOST_TEST(averager.process(std::move(buffer0)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));
  BOOST_TEST(averager.process(std::move(buffer1)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));
  BOOST_TEST(averager.process(std::move(buffer2)));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  Finish(averager, *mock_step);

  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), std::size_t(2));
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[0]->GetRows().size());
  BOOST_REQUIRE_EQUAL(1u, mock_step->GetBdaBuffers()[1]->GetRows().size());
  CheckRow(*average01, *mock_step->GetBdaBuffers()[0], 0, 0);
  CheckRow(*average2, *mock_step->GetBdaBuffers()[1], 0, 0);
}

BOOST_AUTO_TEST_CASE(three_baselines_time_averaging) {
  // This test uses three baselines with lengths of 2400, 300 and 200.
  // The averager should copy the contents of the first, long, baseline.
  // The averager should average the second baseline by a factor of 2.
  // The averager should average the third baseline by a factor of 3.
  const std::size_t kNBaselines = 3;

  const DPInfo info = InitInfo(kAnt1_3Bl, kAnt2_3Bl);
  const double time_threshold =
      600.0;  // Time averaging factors become 1, 2, 3.
  const double chan_threshold = 100.0;  // No channel averaging.

  const dp3::common::ParameterSet parset =
      GetParset(time_threshold, chan_threshold);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  // Create input buffers for the averager. For the first baseline, these
  // buffers also equal the expected output data (expected_long_baseline).
  std::vector<std::unique_ptr<DPBuffer>> buffers;
  std::vector<std::unique_ptr<DPBuffer>> expected_long_baseline;
  for (int i = 0; i < 5; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, kChannelCounts, i * 1000.0));
    expected_long_baseline.push_back(
        std::make_unique<DPBuffer>(*buffers.back()));
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

  // The averager should output BdaBuffers with 2 rows, since the BDA buffers
  // then cover about the same interval (on average) as the input buffers.

  // First baseline: 1 row, second: 0 rows, third: 0 rows -> 0.5 BDA buffers.
  BOOST_TEST(averager.process(std::move(buffers[0])));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(0));

  // First baseline: 2 rows, second: 1 row, third: 0 rows -> 1.5 BDA buffers.
  BOOST_TEST(averager.process(std::move(buffers[1])));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(1));

  // First baseline: 3 rows, second: 1 row, third: 1 row -> 2.5 BDA buffers.
  BOOST_TEST(averager.process(std::move(buffers[2])));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(2));

  // First baseline: 4 rows, second: 2 rows, third: 1 row -> 3.5 BDA buffers.
  BOOST_TEST(averager.process(std::move(buffers[3])));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(3));

  // First baseline: 5 rows, second: 2 rows, third: 1 row -> 4 BDA buffers.
  BOOST_TEST(averager.process(std::move(buffers[4])));
  BOOST_TEST(mock_step->GetBdaBuffers().size() == std::size_t(4));

  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  // First baseline: 5 rows, second: 3 rows, third: 2 rows -> 5 BDA buffers.
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), std::size_t(5));
  for (int i = 0; i < 5; ++i) {
    BOOST_REQUIRE_EQUAL(2u, bda_buffers[i]->GetRows().size());
  }
  CheckRow(*expected_long_baseline[0], *bda_buffers[0], 0, 0);
  CheckRow(*expected_long_baseline[1], *bda_buffers[0], 1, 0);
  CheckRow(*expected_second[0], *bda_buffers[1], 0, 1);
  CheckRow(*expected_long_baseline[2], *bda_buffers[1], 1, 0);
  CheckRow(*expected_third[0], *bda_buffers[2], 0, 2);
  CheckRow(*expected_long_baseline[3], *bda_buffers[2], 1, 0);
  CheckRow(*expected_second[1], *bda_buffers[3], 0, 1);
  CheckRow(*expected_long_baseline[4], *bda_buffers[3], 1, 0);
  CheckRow(*expected_second[2], *bda_buffers[4], 0, 1);
  CheckRow(*expected_third[1], *bda_buffers[4], 1, 2);
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

  const DPInfo info =
      InitInfo(kAnt1_3Bl, kAnt2_3Bl, kInputChannelCounts.size());
  const double time_threshold = 100.0;  // No averaging, baselines are longer.
  const double chan_threshold = 600.0;  // Averaging factors become 1, 2, 3.

  const dp3::common::ParameterSet parset =
      GetParset(time_threshold, chan_threshold);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  // Create input buffers and expected output data for the averager.
  // For the first baseline, the expected output data (expected1) is equal to
  // the input buffer.
  std::vector<std::unique_ptr<DPBuffer>> buffers;
  std::vector<std::unique_ptr<DPBuffer>> expected1;
  std::vector<std::unique_ptr<DPBuffer>> expected2;
  std::vector<std::unique_ptr<DPBuffer>> expected3;
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    const double time = kStartTime + i * kInterval;
    buffers.push_back(CreateBuffer(time, kInterval, kNBaselines,
                                   kInputChannelCounts, i * 1000.0));
    expected1.push_back(std::make_unique<DPBuffer>(*buffers.back()));
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
    BOOST_TEST(averager.process(std::move(buffers[i])));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == i + 1);
  }

  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(0u == bda_buffers[i]->GetRemainingCapacity());
    BOOST_REQUIRE_EQUAL(3u, bda_buffers[i]->GetRows().size());
    CheckRow(*expected1[i], *bda_buffers[i], 0, 0);
    CheckRow(*expected2[i], *bda_buffers[i], 1, 1);
    CheckRow(*expected3[i], *bda_buffers[i], 2, 2);
  }
}

BOOST_AUTO_TEST_CASE(shape_mismatch) {
  // The antenna vectors indicate there is a single baseline.
  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const dp3::common::ParameterSet parset = GetParset(2.0);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  // This buffer has three baselines.
  std::unique_ptr<DPBuffer> buffer =
      CreateBuffer(kStartTime, kInterval, 3, kChannelCounts, 0.0);
  BOOST_CHECK_THROW(averager.process(std::move(buffer)), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(max_interval) {
  const double kMaxInterval = kInterval * 3.5;
  const std::size_t kFactor = std::floor(kMaxInterval);
  const std::size_t kNBaselines = 1;
  const std::size_t kTimeSteps = 7 * kFactor;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  // The threshold implies an averaging factor of 42, however, kMaxInterval
  // should limit the averaging to a factor of 3.
  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * 42.0, std::nullopt, kMaxInterval);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(
        averager.process(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                      kNBaselines, kChannelCounts, 0.0)));
    BOOST_TEST(mock_step->GetBdaBuffers().size() == (i + 1) / kFactor);
  }

  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), kTimeSteps / kFactor);
  for (const auto& bda_buffer : bda_buffers) {
    const std::vector<BdaBuffer::Row> rows = bda_buffer->GetRows();
    BOOST_REQUIRE_EQUAL(1u, rows.size());
    BOOST_TEST(kInterval * kFactor == rows[0].interval);
  }
}

BOOST_AUTO_TEST_CASE(min_channels) {
  const std::size_t kMinChannels = 3;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kOutputChannelCounts{1, 2, 2};

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  // The threshold implies that all channels should be averaged into one
  // channel, however, we specify a minimum of three channels per baseline.
  const dp3::common::ParameterSet parset = GetParset(
      std::nullopt, baseline_length * kNChan, std::nullopt, kMinChannels);
  BdaAverager averager(parset, kPrefix);
  BOOST_REQUIRE_NO_THROW(averager.updateInfo(info));

  std::unique_ptr<DPBuffer> buffer =
      CreateBuffer(kStartTime, kInterval, kNBaselines, kChannelCounts, 0.0);
  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(std::move(buffer)));
  Finish(averager, *mock_step);

  const auto& bda_buffers = mock_step->GetBdaBuffers();
  BOOST_REQUIRE_EQUAL(bda_buffers.size(), std::size_t(1));
  CheckRow(*averaged, *bda_buffers[0], 0, 0);
}

BOOST_AUTO_TEST_CASE(zero_values_weight) {
  const std::size_t kMinChannels = 3;
  const std::size_t kNBaselines = 1;
  const std::vector<std::size_t> kOutputChannelCounts{1, 2, 2};

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  // The threshold implies that all channels should be averaged into one
  // channel, however, we specify a minimum of three channels per baseline.
  const dp3::common::ParameterSet parset = GetParset(
      std::nullopt, baseline_length * kNChan, std::nullopt, kMinChannels);
  BdaAverager averager(parset, kPrefix);
  averager.updateInfo(info);

  // Calling CreateBuffer with a weight of 0.0 yields NaN data values.
  std::unique_ptr<DPBuffer> buffer = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kChannelCounts, 0.0, 0.0);
  // Put non-NaN values into the imaginary part of the data.
  for (std::complex<float>& value : buffer->GetData()) {
    value = {value.real(), 42.0f};
  }
  // Put non-NaN values into the real part of the extra data.
  for (std::complex<float>& value : buffer->GetData(kExtraData)) {
    value = {0.42f, value.imag()};
  }

  std::unique_ptr<DPBuffer> averaged = CreateBuffer(
      kStartTime, kInterval, kNBaselines, kOutputChannelCounts, 0.0);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  averager.setNextStep(mock_step);

  BOOST_TEST(averager.process(std::move(buffer)));
  Finish(averager, *mock_step);

  // A weight of 0 should only result in NaN when the input data is NaN.
  // When the input data is not NaN, the result should be zero.
  const auto& bda_buffers = mock_step->GetBdaBuffers();

  const size_t data_size = bda_buffers[0]->GetNumberOfElements();
  BOOST_TEST(data_size > 0);
  for (std::size_t i = 0; i < data_size; ++i) {
    const std::complex<float> value = bda_buffers[0]->GetData()[i];
    const std::complex<float> extra_value =
        bda_buffers[0]->GetData(kExtraData)[i];
    BOOST_TEST(std::isnan(value.real()));
    BOOST_TEST(value.imag() == 0);
    BOOST_TEST(extra_value.real() == 0);
    BOOST_TEST(std::isnan(extra_value.imag()));
  }
}

BOOST_AUTO_TEST_CASE(force_buffersize) {
  const std::size_t kTimeAveragingFactor = 2;
  const std::size_t kNBaselines = 1;
  const int kNInputBuffers = 20;

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * kTimeAveragingFactor);
  BdaAverager averager(parset, kPrefix);
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
    BOOST_TEST(averager.process(std::move(input_buffers[i])));
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

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * kTimeAveragingFactor);
  BdaAverager averager(parset, kPrefix);
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
    BOOST_TEST(averager.process(std::move(input_buffers[i])));
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

  const DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  const double baseline_length = info.getBaselineLengths()[0];

  const dp3::common::ParameterSet parset =
      GetParset(baseline_length * kTimeAveragingFactor);
  BdaAverager averager(parset, kPrefix);
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
    BOOST_TEST(averager.process(std::move(input_buffers[i])));
  }

  Finish(averager, *mock_step);

  // With the given input, we expect a total of 200 elements in the output. In
  // this test we check that output has 200 elements in case we force
  // the output buffersize to a value higher than the maximum (all the elements
  // are forced into one single BdaBuffer).
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers().size(), 1u);
  BOOST_REQUIRE_EQUAL(mock_step->GetBdaBuffers()[0]->GetNumberOfElements(),
                      200u);
}

BOOST_AUTO_TEST_SUITE_END()
