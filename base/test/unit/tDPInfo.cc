// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../DPInfo.h"

#include <boost/test/unit_test.hpp>

namespace {
const std::vector<std::string> kAntNames{"foobar"};
const std::vector<double> kAntDiam = {42.0};
const std::vector<casacore::MPosition> kAntPos(1);
}  // namespace

BOOST_AUTO_TEST_SUITE(dpinfo)

BOOST_AUTO_TEST_CASE(set_frequency_info) {
  const std::vector<double> kFreqs{10.0, 20.0};
  const std::vector<double> kWidths{5.0, 6.0};
  const double kRefFreq = 15.0;
  const double kTotalWidth = 11.0;

  dp3::base::DPInfo info;
  info.set(std::vector<double>(kFreqs), std::vector<double>(kWidths));
  BOOST_TEST(kFreqs == info.chanFreqs());
  BOOST_TEST(kWidths == info.chanWidths());
  BOOST_TEST(kWidths == info.resolutions());
  BOOST_TEST(kWidths == info.effectiveBW());
  BOOST_TEST(kRefFreq == info.refFreq());
  BOOST_TEST(kTotalWidth == info.totalBW());
}

BOOST_AUTO_TEST_CASE(set_bda_frequency_info) {
  const std::vector<int> kAnt(3, 0);
  const std::vector<std::vector<double>> kFreqs{
      {30.0}, {10.0, 20.0, 30.0, 40.0, 50.0}, {20.0, 45.0}};
  const std::vector<std::vector<double>> kWidths{
      {50.0}, {10.0, 10.0, 10.0, 10.0, 10.0}, {30.0, 20.0}};
  const double kRefFreq = 30.0;
  const double kTotalWidth = 50.0;

  dp3::base::DPInfo info;
  info.set(kAntNames, kAntDiam, kAntPos, kAnt, kAnt);  // Set baseline count.
  info.set(std::vector<std::vector<double>>(kFreqs),
           std::vector<std::vector<double>>(kWidths));
  for (std::size_t i = 0; i < kFreqs.size(); i++) {
    BOOST_TEST(kFreqs[i] == info.chanFreqs(i));
    BOOST_TEST(kWidths[i] == info.chanWidths(i));
    BOOST_TEST(kWidths[i] == info.resolutions(i));
    BOOST_TEST(kWidths[i] == info.effectiveBW(i));
  }
  BOOST_TEST(kRefFreq == info.refFreq());
  BOOST_TEST(kTotalWidth == info.totalBW());
}

BOOST_AUTO_TEST_CASE(channels_are_regular) {
  // Note that the tolerance in channelsAreRegular is 1000 Hz.

  const std::vector<int> kAnt(3, 0);

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
    dp3::base::DPInfo info;
    info.set(std::vector<double>(kRegularFreqs.front()),
             std::vector<double>(kRegularWidths.front()));
    BOOST_TEST(info.channelsAreRegular());
  }
  {
    dp3::base::DPInfo info;
    info.set(std::vector<double>(kRegularFreqs.front()),
             std::vector<double>(kIrregularWidths));
    BOOST_TEST(!info.channelsAreRegular());
  }
  {
    dp3::base::DPInfo info;
    info.set(std::vector<double>(kIrregularFreqs),
             std::vector<double>(kRegularWidths.front()));
    BOOST_TEST(!info.channelsAreRegular());
  }

  // Test using multiple baselines.
  {
    dp3::base::DPInfo info;
    info.set(kAntNames, kAntDiam, kAntPos, kAnt, kAnt);  // Set baseline count.
    info.set(std::vector<std::vector<double>>(kRegularFreqs),
             std::vector<std::vector<double>>(kRegularWidths));
    BOOST_TEST(info.channelsAreRegular());
  }
  {
    dp3::base::DPInfo info;
    info.set(kAntNames, kAntDiam, kAntPos, kAnt, kAnt);  // Set baseline count.
    info.set(std::vector<std::vector<double>>(kIrregularFreqsBDA),
             std::vector<std::vector<double>>(kIrregularWidthsBDA));
    BOOST_TEST(!info.channelsAreRegular());
  }
}

BOOST_AUTO_TEST_SUITE_END()
