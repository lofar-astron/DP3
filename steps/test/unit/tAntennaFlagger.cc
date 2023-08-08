// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <random>
#include <xtensor/xview.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xrandom.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include "../../AntennaFlagger.h"
#include "../../../common/baseline_indices.h"
#include "../../../common/ParameterSet.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::AntennaFlagger;

namespace {

static const size_t kBrokenAntenna = 42;
static const size_t kBrokenStation = 6;
static const float kFirstTime = 100.0;
static const float kTimeStep = 5.0;
static const float kStartFrequency = 1050000.0;
static const float kChannelWidth = 1.0;

/**
 * Class to generate random test input data with some simulated broken antennas:
 * - For one station, one antenna is broken. Set the data for all affected
 *   baselines to some arbitrary value far outside of the range of the random
 *   input data to make sure it stands out.
 * - For another station, all antennas are broken. 'Corrupt' the data for the
 *   affected baselines a bit by multiplying the data with a factor of 2x.
 *
 * The goal of this test is to have artificial data in which both the single
 * broken antenna, as well as all antennas in the broken station can be detected
 * by the flagger.
 */
class TestInput final : public dp3::steps::MockInput {
 public:
  TestInput(size_t n_time, size_t n_stations, size_t n_antennas_per_station,
            size_t n_channels, size_t n_correlations,
            const std::string& antenna_name_prefix,
            dp3::common::BaselineOrder baseline_order)
      : n_time_processed_(0),
        n_time_(n_time),
        n_stations_(n_stations),
        n_antennas_per_station_(n_antennas_per_station),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        baseline_order_(baseline_order) {
    info() = DPInfo(n_correlations, n_channels);
    info().setTimes(kFirstTime, kFirstTime + (n_time_ - 1) * kTimeStep,
                    kTimeStep);

    // Fill the baselines
    const size_t n_antennas = n_stations * n_antennas_per_station;
    n_baselines_ = dp3::common::ComputeNBaselines(n_antennas);
    std::vector<int> ant1(n_baselines_);
    std::vector<int> ant2(n_baselines_);

    for (size_t i = 0; i < n_baselines_; ++i) {
      const std::pair<size_t, size_t> antennas =
          dp3::common::ComputeBaseline(i, n_antennas, baseline_order);
      ant1[i] = antennas.first;
      ant2[i] = antennas.second;
    }

    std::vector<std::string> antenna_names(n_antennas);
    for (size_t i = 0; i < n_antennas; ++i) {
      antenna_names[i] = antenna_name_prefix + std::to_string(i);
    }

    // Antenna positions and diameter are not used by the AntennaFlagger
    std::vector<casacore::MPosition> antenna_positions(n_antennas);
    std::vector<double> antenna_diameters(n_antennas);
    info().setAntennas(antenna_names, antenna_diameters, antenna_positions,
                       ant1, ant2);

    // Define the frequencies
    std::vector<double> frequencies(n_channels_);
    std::vector<double> channel_widths(n_channels_, kChannelWidth);
    std::iota(frequencies.begin(), frequencies.end(), kStartFrequency);
    info().setChannels(std::move(frequencies), std::move(channel_widths));
  }

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (n_time_processed_ == n_time_) {
      return false;
    }

    // Initialize data
    buffer->GetData().resize({n_baselines_, n_channels_, n_correlations_});
    xt::random::seed(0);
    xt::real(buffer->GetData()) =
        xt::random::randn<float>(buffer->GetData().shape(), 0, 1);
    xt::imag(buffer->GetData()) =
        xt::random::randn<float>(buffer->GetData().shape(), 0, 1);

    // Set broken antenna
    const size_t n_antennas = n_stations_ * n_antennas_per_station_;

    for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
      const size_t bl = dp3::common::ComputeBaselineIndex(
          kBrokenAntenna, antenna2, n_antennas, baseline_order_);
      xt::view(buffer->GetData(), bl, xt::all(), xt::all()).fill(42.0f);
    }

    // Set broken station
    for (size_t i = 0; i < n_antennas_per_station_; ++i) {
      const size_t antenna1 = kBrokenStation * n_antennas_per_station_ + i;
      assert(antenna1 != kBrokenAntenna);
      for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
        const size_t bl = dp3::common::ComputeBaselineIndex(
            antenna1, antenna2, n_antennas, baseline_order_);
        xt::view(buffer->GetData(), bl, xt::all(), xt::all()) *= 2.0f;
      }
    }

    getNextStep()->process(std::move(buffer));
    ++n_time_processed_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void show(std::ostream&) const override{};
  void updateInfo(const DPInfo&) override {}

  size_t n_time_processed_;
  size_t n_time_;
  size_t n_stations_;
  size_t n_antennas_per_station_;
  size_t n_baselines_;
  size_t n_channels_;
  size_t n_correlations_;
  dp3::common::BaselineOrder baseline_order_;
};

/**
 * Class to verify that the 'broken' antennas are indeed flagged. No other
 * antennas should be flagged.
 */
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(size_t n_time, size_t n_baselines, size_t n_stations,
             size_t n_antennas_per_station, size_t n_channels,
             size_t n_correlations, dp3::common::BaselineOrder baseline_order)
      : n_time_(n_time),
        n_baselines_(n_baselines),
        n_stations_(n_stations),
        n_antennas_per_station_(n_antennas_per_station),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        baseline_order_(baseline_order) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    xt::xtensor<bool, 1> result({n_baselines_}, false);

    const size_t n_antennas = n_stations_ * n_antennas_per_station_;
    for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
      const size_t bl = dp3::common::ComputeBaselineIndex(
          kBrokenAntenna, antenna2, n_antennas, baseline_order_);
      result(bl) = true;
    }

    for (size_t i = 0; i < n_antennas_per_station_; ++i) {
      const size_t antenna1 = kBrokenStation * n_antennas_per_station_ + i;
      assert(antenna1 != kBrokenAntenna);
      for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
        const size_t bl = dp3::common::ComputeBaselineIndex(
            antenna1, antenna2, n_antennas, baseline_order_);
        result(bl) = true;
      }
    }

    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      if (result(bl)) {
        BOOST_CHECK(
            xt::all(xt::view(buffer->GetFlags(), bl, xt::all(), xt::all())));
      } else {
        BOOST_CHECK(
            !xt::any(xt::view(buffer->GetFlags(), bl, xt::all(), xt::all())));
      }
    }

    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(static_cast<int>(info.origNChan()), n_channels_);
    BOOST_CHECK_EQUAL(static_cast<int>(info.nchan()), n_channels_);
    BOOST_CHECK(info.ntime() == n_time_);
    BOOST_CHECK(info.firstTime() == kFirstTime);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5);
    BOOST_CHECK_EQUAL(static_cast<int>(info.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(static_cast<int>(info.ntimeAvg()), 1);
    BOOST_CHECK(info.chanFreqs().size() == n_channels_);
    BOOST_CHECK(info.chanWidths().size() == n_channels_);
  }

  size_t n_time_;
  size_t n_baselines_;
  size_t n_stations_;
  size_t n_antennas_per_station_;
  size_t n_channels_;
  size_t n_correlations_;
  dp3::common::BaselineOrder baseline_order_;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(antennaflagger)

BOOST_AUTO_TEST_CASE(aartfaac) {
  const int kNTimes = 1;
  const int kNStations = 12;
  const int kNAntennasPerStation = 48;
  const int kNChannels = 8;
  const int kNCorrelations = 4;
  const std::string kAntennaNamePrefix = "A12_";
  const dp3::common::BaselineOrder kBaselineOrder =
      dp3::common::BaselineOrder::kRowMajor;
  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  auto in = std::make_shared<TestInput>(
      kNTimes, kNStations, kNAntennasPerStation, kNChannels, kNCorrelations,
      kAntennaNamePrefix, kBaselineOrder);
  dp3::common::ParameterSet parset;
  parset.add("antenna_flagging_sigma", "4");
  parset.add("antenna_flagging_max_iterations", "1");
  parset.add("station_flagging_sigma", "3");
  parset.add("station_flagging_max_iterations", "1");
  auto antenna_flagger = std::make_shared<AntennaFlagger>(parset, "");
  BOOST_CHECK(parset.unusedKeys().empty());
  auto out = std::make_shared<TestOutput>(kNTimes, kNBaselines, kNStations,
                                          kNAntennasPerStation, kNChannels,
                                          kNCorrelations, kBaselineOrder);
  dp3::steps::test::Execute({in, antenna_flagger, out});
}

BOOST_AUTO_TEST_CASE(lofar) {
  const int kNTimes = 1;
  const int kNStations = 55;
  const int kNAntennasPerStation = 1;
  const int kNChannels = 8;
  const int kNCorrelations = 4;
  const std::string kAntennaNamePrefix = "";
  const dp3::common::BaselineOrder kBaselineOrder =
      dp3::common::BaselineOrder::kColumnMajor;
  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  auto in = std::make_shared<TestInput>(
      kNTimes, kNStations, kNAntennasPerStation, kNChannels, kNCorrelations,
      kAntennaNamePrefix, kBaselineOrder);
  dp3::common::ParameterSet parset;
  parset.add("antenna_flagging_sigma", "3");
  parset.add("antenna_flagging_max_iterations", "1");
  parset.add("station_flagging_sigma", "2");
  parset.add("station_flagging_max_iterations", "1");
  auto antenna_flagger = std::make_shared<AntennaFlagger>(parset, "");
  BOOST_CHECK(parset.unusedKeys().empty());
  auto out = std::make_shared<TestOutput>(kNTimes, kNBaselines, kNStations,
                                          kNAntennasPerStation, kNChannels,
                                          kNCorrelations, kBaselineOrder);
  dp3::steps::test::Execute({in, antenna_flagger, out});
}

BOOST_AUTO_TEST_CASE(incorrect_baseline_count) {
  const int kNStations = 2;
  const int kNAntennasPerStation = 2;
  const int kNChannels = 8;
  const int kNCorrelations = 4;
  const int kNAntennas = kNStations * kNAntennasPerStation;

  // Initialize DPInfo with only two baselines, while the AntennaFlagger expects
  // a baseline for every antenna pair.
  const std::vector<int> ant1{0, 0};
  const std::vector<int> ant2{0, 1};

  std::vector<std::string> antenna_names(kNAntennas);
  for (size_t i = 0; i < kNAntennas; ++i) {
    antenna_names[i] = std::to_string(i);
  }

  DPInfo info(kNCorrelations, kNChannels);
  const std::vector<casacore::MPosition> antenna_positions(kNAntennas);
  const std::vector<double> antenna_diameters(kNAntennas);
  info.setAntennas(antenna_names, antenna_diameters, antenna_positions, ant1,
                   ant2);

  // Initialize AntennaFlagger
  dp3::common::ParameterSet parset;
  auto antenna_flagger = std::make_shared<AntennaFlagger>(parset, "");

  // Calling ::updateInfo should throw an exception
  BOOST_CHECK_THROW(antenna_flagger->updateInfo(info), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(incorrect_baseline_order) {
  const int kNStations = 2;
  const int kNAntennasPerStation = 2;
  const int kNChannels = 8;
  const int kNCorrelations = 4;
  const int kNAntennas = kNStations * kNAntennasPerStation;
  const int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);
  const dp3::common::BaselineOrder kBaselineOrder =
      dp3::common::BaselineOrder::kRowMajor;

  // Initialize DPInfo with antennas swapped
  std::vector<int> ant1(kNBaselines);
  std::vector<int> ant2(kNBaselines);

  for (int i = 0; i < kNBaselines; ++i) {
    const std::pair<size_t, size_t> antennas =
        dp3::common::ComputeBaseline(i, kNAntennas, kBaselineOrder);
    ant1[i] = antennas.second;  // should be antennas.first
    ant2[i] = antennas.first;   // should be antennas.second
  }

  std::vector<std::string> antenna_names(kNAntennas);
  for (int i = 0; i < kNAntennas; ++i) {
    antenna_names[i] = std::to_string(i);
  }

  DPInfo info(kNCorrelations, kNChannels);
  const std::vector<casacore::MPosition> antenna_positions(kNAntennas);
  const std::vector<double> antenna_diameters(kNAntennas);
  info.setAntennas(antenna_names, antenna_diameters, antenna_positions, ant1,
                   ant2);

  // Initialize AntennaFlagger
  dp3::common::ParameterSet parset;
  auto antenna_flagger = std::make_shared<AntennaFlagger>(parset, "");

  // Calling ::updateInfo should throw an exception
  BOOST_CHECK_THROW(antenna_flagger->updateInfo(info), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
