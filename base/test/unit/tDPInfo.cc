// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "base/DPInfo.h"
#include <boost/test/unit_test.hpp>

using dp3::base::DPInfo;

namespace {
const std::size_t kNBaselines = 3;
const std::vector<int> kAntennaIndices(kNBaselines, 0);
const std::vector<std::string> kAntennaNames{"foobar"};
const std::vector<double> kAntennaDiameters = {42.0};
const std::vector<casacore::MPosition> kAntennaPositions(1);
}  // namespace

BOOST_AUTO_TEST_SUITE(dpinfo)

BOOST_AUTO_TEST_CASE(constructor) {
  const DPInfo default_values;
  BOOST_TEST(default_values.ncorr() == 0);
  BOOST_TEST(default_values.origNChan() == 0);
  BOOST_TEST(default_values.nchan() == 0);
  BOOST_TEST(default_values.antennaSet().empty());

  const unsigned int kNCorrelations = 4;
  const unsigned int kNChannels = 42;
  const std::string kAntennaSet = "test_antenna_set";
  const DPInfo custom_values(kNCorrelations, kNChannels, kAntennaSet);
  BOOST_TEST(custom_values.ncorr() == kNCorrelations);
  BOOST_TEST(custom_values.origNChan() == kNChannels);
  BOOST_TEST(custom_values.nchan() == kNChannels);
  BOOST_TEST(custom_values.antennaSet() == kAntennaSet);
}

BOOST_AUTO_TEST_CASE(set_times) {
  DPInfo info;
  // Check default values.
  BOOST_TEST(info.firstTime() == 0.0);
  BOOST_TEST(info.lastTime() == 0.0);
  BOOST_TEST(info.timeInterval() == 1.0);
  BOOST_TEST(info.startTime() == -0.5);
  BOOST_TEST(info.ntime() == 1);

  const double kFirstTime = 42.0;
  const double kLastTime = 49.0;
  const double kInterval = 4.0;
  info.setTimes(kFirstTime, kLastTime, kInterval);
  BOOST_TEST(info.firstTime() == kFirstTime);
  BOOST_TEST(info.lastTime() == kLastTime);
  BOOST_TEST(info.timeInterval() == kInterval);
  BOOST_TEST(info.startTime() == 40.0);
  BOOST_TEST(info.ntime() == 3);
}

BOOST_AUTO_TEST_CASE(set_negative_times) {
  const double kFirstTime = -100.0;
  const double kLastTime = -50.0;
  const double kInterval = 2.0;
  DPInfo info;
  info.setTimes(kFirstTime, kLastTime, kInterval);
  BOOST_TEST(info.firstTime() == kFirstTime);
  BOOST_TEST(info.lastTime() == kLastTime);
  BOOST_TEST(info.timeInterval() == kInterval);
  BOOST_TEST(info.startTime() == -101.0);
  BOOST_TEST(info.ntime() == 26);
}

BOOST_AUTO_TEST_CASE(set_times_last_before_first) {
  const double kFirstTime = 42.0;
  const double kLastTime = 0.0;
  const double kInterval = 1.0;
  DPInfo info;
  BOOST_CHECK_THROW(info.setTimes(kFirstTime, kLastTime, kInterval),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_times_invalid_interval) {
  const double kFirstTime = 0.0;
  const double kLastTime = 1.0;
  const double kZeroInterval = 0.0;
  const double kNegativeInterval = -42.0;
  {
    DPInfo info;
    BOOST_CHECK_THROW(info.setTimes(kFirstTime, kLastTime, kZeroInterval),
                      std::invalid_argument);
  }
  {
    DPInfo info;
    BOOST_CHECK_THROW(info.setTimes(kFirstTime, kLastTime, kNegativeInterval),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(set_array_info) {
  const casacore::Vector<double> kPositionValues{
      std::vector<int>{4'200'000, 420'000, 4'200'042}};
  const casacore::MPosition kArrayPosition(
      casacore::Quantum(kPositionValues, "m"), casacore::MPosition::ITRF);
  const casacore::MDirection kPhaseCenter(casacore::Quantity(42, "deg"),
                                          casacore::Quantity(1, "deg"),
                                          casacore::MDirection::J2000);
  const casacore::MDirection kDelayCenter(casacore::Quantity(43, "deg"),
                                          casacore::Quantity(2, "deg"),
                                          casacore::MDirection::J2000);
  const casacore::MDirection kTimeBeamDirection(casacore::Quantity(44, "deg"),
                                                casacore::Quantity(3, "deg"),
                                                casacore::MDirection::J2000);

  DPInfo info;
  BOOST_TEST(info.arrayPos().getValue() == casacore::MPosition().getValue());
  BOOST_TEST(info.originalPhaseCenter().getValue() ==
             casacore::MDirection().getValue());
  BOOST_TEST(info.phaseCenter().getValue() ==
             casacore::MDirection().getValue());
  BOOST_TEST(info.delayCenter().getValue() ==
             casacore::MDirection().getValue());
  BOOST_TEST(info.tileBeamDir().getValue() ==
             casacore::MDirection().getValue());

  info.setArrayInformation(kArrayPosition, kPhaseCenter, kDelayCenter,
                           kTimeBeamDirection);
  BOOST_TEST(info.arrayPos().getValue() == kArrayPosition.getValue());
  BOOST_TEST(info.originalPhaseCenter().getValue() == kPhaseCenter.getValue());
  BOOST_TEST(info.phaseCenter().getValue() == kPhaseCenter.getValue());
  BOOST_TEST(info.delayCenter().getValue() == kDelayCenter.getValue());
  BOOST_TEST(info.tileBeamDir().getValue() == kTimeBeamDirection.getValue());
}

BOOST_AUTO_TEST_CASE(set_channels) {
  const std::vector<double> kFreqs{10.0, 20.0};
  const std::vector<double> kWidths{5.0, 6.0};
  const double kRefFreq = 15.0;
  const double kTotalWidth = 11.0;
  const unsigned int kNOriginalChannels = 42;

  DPInfo info(1, kNOriginalChannels);
  BOOST_TEST(info.origNChan() == kNOriginalChannels);
  info.setChannels(std::vector<double>(kFreqs), std::vector<double>(kWidths));
  // setChannels() should set change the original channel count.
  BOOST_TEST(kNOriginalChannels == info.origNChan());
  BOOST_TEST(kFreqs.size() == info.nchan());
  BOOST_TEST(kFreqs == info.chanFreqs());
  BOOST_TEST(kWidths == info.chanWidths());
  BOOST_TEST(kWidths == info.resolutions());
  BOOST_TEST(kWidths == info.effectiveBW());
  BOOST_TEST(kRefFreq == info.refFreq());
  BOOST_TEST(kTotalWidth == info.totalBW());
}

BOOST_AUTO_TEST_CASE(select_channels) {
  DPInfo info;
  BOOST_CHECK(info.startchan() == 0);
  BOOST_CHECK(info.nchan() == 0);

  // Create two channels.
  const std::vector<double> kFreqs{10.0, 20.0};
  const std::vector<double> kWidths{5.0, 6.0};
  info.setChannels(std::vector<double>(kFreqs), std::vector<double>(kWidths));
  BOOST_CHECK(info.startchan() == 0);
  BOOST_CHECK(info.nchan() == kFreqs.size());

  // Select only the second channel.
  const unsigned int kStartChannel = 1;
  const unsigned int kNChannel = 1;
  BOOST_CHECK_NO_THROW(info.SelectChannels(kStartChannel, kNChannel));
  BOOST_CHECK(info.startchan() == kStartChannel);
  BOOST_CHECK(info.nchan() == kNChannel);

  // Try to select only the second channel again.
  // This should fail because only one channel is left.
  BOOST_CHECK_THROW(info.SelectChannels(kStartChannel, kNChannel),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_bda_channels) {
  const std::vector<std::vector<double>> kFreqs{
      {30.0}, {10.0, 20.0, 30.0, 40.0, 50.0}, {20.0, 45.0}};
  const std::vector<std::vector<double>> kWidths{
      {50.0}, {10.0, 10.0, 10.0, 10.0, 10.0}, {30.0, 20.0}};
  const double kRefFreq = 30.0;
  const double kTotalWidth = 50.0;
  const unsigned int kNOriginalChannels = 42;

  DPInfo info(1, kNOriginalChannels);
  info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                   kAntennaIndices,
                   kAntennaIndices);  // Set baseline count.
  info.setChannels(std::vector<std::vector<double>>(kFreqs),
                   std::vector<std::vector<double>>(kWidths));
  for (std::size_t i = 0; i < kFreqs.size(); i++) {
    BOOST_TEST(kFreqs[i] == info.chanFreqs(i));
    BOOST_TEST(kWidths[i] == info.chanWidths(i));
    BOOST_TEST(kWidths[i] == info.resolutions(i));
    BOOST_TEST(kWidths[i] == info.effectiveBW(i));
  }
  // setChannels() should set change the original channel count.
  BOOST_TEST(kNOriginalChannels == info.origNChan());
  BOOST_TEST(5 == info.nchan());  // Maximum number of channels.
  BOOST_TEST(kRefFreq == info.refFreq());
  BOOST_TEST(kTotalWidth == info.totalBW());
}

BOOST_AUTO_TEST_CASE(set_mismatching_channels) {
  std::vector<double> kFrequencies{42.0};
  std::vector<double> kWidths{1.0};
  std::vector<double> kMismatch{1.0, 2.0};

  BOOST_CHECK_NO_THROW(DPInfo().setChannels(std::vector<double>(kFrequencies),
                                            std::vector<double>(kWidths)));

  // Use kMismatch for widths.
  BOOST_CHECK_THROW(DPInfo().setChannels(std::vector<double>(kFrequencies),
                                         std::vector<double>(kMismatch)),
                    std::invalid_argument);

  // Use kMisMatch for resolutions.
  BOOST_CHECK_THROW(DPInfo().setChannels(std::vector<double>(kFrequencies),
                                         std::vector<double>(kWidths),
                                         std::vector<double>(kMismatch)),
                    std::invalid_argument);

  // Use kMisMatch for effective bandwidth.
  BOOST_CHECK_THROW(DPInfo().setChannels(std::vector<double>(kFrequencies),
                                         std::vector<double>(kWidths), {},
                                         std::vector<double>(kMismatch)),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_mismatching_bda_channels) {
  const std::vector<std::vector<double>> kCorrectFrequencies{
      {42.0}, {42.0}, {42.0}};
  const std::vector<std::vector<double>> kCorrectWidths{{1.0}, {1.0}, {1.0}};

  DPInfo info;
  info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                   kAntennaIndices,
                   kAntennaIndices);  // Set baseline count.
  {
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths = kCorrectWidths;
    BOOST_CHECK_NO_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths)));
  }

  {  // Use an incorrect baseline count for frequencies.
    std::vector<std::vector<double>> incorrect_frequencies{{42.0}, {42.0}};
    std::vector<std::vector<double>> widths = kCorrectWidths;
    BOOST_CHECK_THROW(DPInfo(info).setChannels(std::move(incorrect_frequencies),
                                               std::move(widths)),
                      std::invalid_argument);
  }
  {  // Use an incorrect baseline count for widths.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> incorrect_widths{{1.0}, {1.0}};
    BOOST_CHECK_THROW(DPInfo(info).setChannels(std::move(frequencies),
                                               std::move(incorrect_widths)),
                      std::invalid_argument);
  }
  {  // Use an incorrect baseline count for resolutions.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths = kCorrectWidths;
    std::vector<std::vector<double>> incorrect_resolutions{{1.0}, {1.0}};
    BOOST_CHECK_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths),
                                 std::move(incorrect_resolutions)),
        std::invalid_argument);
  }
  {  // Use an incorrect baseline count for effective bandwidths.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths = kCorrectWidths;
    std::vector<std::vector<double>> incorrect_effective_bandwidths{{1.0},
                                                                    {1.0}};
    BOOST_CHECK_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths), {},
                                 std::move(incorrect_effective_bandwidths)),
        std::invalid_argument);
  }
  {  // Use more widths than frequencies in the first baseline.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths{{1.0, 1.0}, {1.0}, {1.0}};
    BOOST_CHECK_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths)),
        std::invalid_argument);
  }
  {  // Use more resolutions than frequencies in the second baseline.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths = kCorrectWidths;
    std::vector<std::vector<double>> resolutions{{1.0}, {1.0, 1.0}, {1.0}};
    BOOST_CHECK_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths),
                                 std::move(resolutions)),
        std::invalid_argument);
  }
  {  // Use more effective bandwidths than frequencies in the third baseline.
    std::vector<std::vector<double>> frequencies = kCorrectFrequencies;
    std::vector<std::vector<double>> widths = kCorrectWidths;
    std::vector<std::vector<double>> effective_bandwidths{
        {1.0}, {1.0}, {1.0, 1.0}};
    BOOST_CHECK_THROW(
        DPInfo(info).setChannels(std::move(frequencies), std::move(widths), {},
                                 std::move(effective_bandwidths)),
        std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(channels_are_regular) {
  // Note that the tolerance in channelsAreRegular is 1000 Hz.

  const std::vector<std::vector<double>> kRegularFreqs(
      3, {10000.0, 20000.0, 30000.0, 40000.0});
  const std::vector<std::vector<double>> kRegularWidths(
      3, {5000.0, 5000.0, 5000.0, 5000.0});

  const std::vector<double> kIrregularFreqs{10000.0, 20000.0, 30000.0, 42000.0};
  const std::vector<double> kIrregularWidths{5000.0, 7500.0, 5000.0, 5000.0};

  // In the BDA data, each baseline is regular, however, BDA makes it irregular.
  const std::vector<std::vector<double>> kIrregularFreqsBDA{
      {10000.0, 20000.0, 30000.0, 40000.0}, {250000.0}, {250000.0}};
  const std::vector<std::vector<double>> kIrregularWidthsBDA{
      {5000.0, 5000.0, 5000.0, 5000.0}, {20000.0}, {20000.0}};

  // Test using a single baseline.
  {
    DPInfo info;
    info.setChannels(std::vector<double>(kRegularFreqs.front()),
                     std::vector<double>(kRegularWidths.front()));
    BOOST_TEST(info.channelsAreRegular());
  }
  {
    DPInfo info;
    info.setChannels(std::vector<double>(kRegularFreqs.front()),
                     std::vector<double>(kIrregularWidths));
    BOOST_TEST(!info.channelsAreRegular());
  }
  {
    DPInfo info;
    info.setChannels(std::vector<double>(kIrregularFreqs),
                     std::vector<double>(kRegularWidths.front()));
    BOOST_TEST(!info.channelsAreRegular());
  }

  // Test using multiple baselines.
  {
    DPInfo info;
    info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                     kAntennaIndices, kAntennaIndices);  // Set baseline count.
    info.setChannels(std::vector<std::vector<double>>(kRegularFreqs),
                     std::vector<std::vector<double>>(kRegularWidths));
    BOOST_TEST(info.channelsAreRegular());
  }
  {
    DPInfo info;
    info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                     kAntennaIndices, kAntennaIndices);  // Set baseline count.
    info.setChannels(std::vector<std::vector<double>>(kIrregularFreqsBDA),
                     std::vector<std::vector<double>>(kIrregularWidthsBDA));
    BOOST_TEST(!info.channelsAreRegular());
  }
}

BOOST_AUTO_TEST_SUITE_END()
