// tAntennaFlagger.cc: Test program for class AntennaFlagger
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include "../../AntennaFlagger.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::AntennaFlagger;

BOOST_AUTO_TEST_SUITE(antennaflagger)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput final : public dp3::steps::MockInput {
 public:
  TestInput(int n_time, int n_baselines, int n_channels, int n_correlations,
            bool flag)
      : n_time_processed_(0),
        n_time_(n_time),
        n_baselines_(n_baselines),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        flag_(flag) {
    info() = DPInfo(n_correlations, n_channels);
    info().setTimes(0.5, 0.5 + (n_time - 1) * 5.0, 5.0);

    // Fill the baselines; use 12 stations with 48 antennas each.
    constexpr int n_receivers = 12 * 48;
    std::vector<int> ant1(n_baselines);
    std::vector<int> ant2(n_baselines);
    int station1 = 0;
    int station2 = 0;
    for (int i = 0; i < n_baselines; ++i) {
      ant1[i] = station1;
      ant2[i] = station2;
      if (++station2 == n_receivers) {
        station2 = 0;
        if (++station1 == n_receivers) {
          station1 = 0;
        }
      }
    }

    // Atenna names and positions are not used by the AntennaFlagger,
    // but the arrays need to have the correct dimensions.
    std::vector<std::string> antNames(n_receivers);
    std::vector<casacore::MPosition> antPos(n_receivers);
    std::vector<double> antDiam(n_receivers, 70.);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
  }

 private:
  bool process(const DPBuffer&) override {
    // Stop when all times are done.
    if (n_time_processed_ == n_time_) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(n_correlations_, n_channels_,
                                           n_baselines_);
    for (size_t i = 0; i < data.size(); ++i) {
      data.data()[i] = casacore::Complex(i + n_time_processed_ * 10,
                                         i - 10 + n_time_processed_ * 6);
    }
    casacore::Matrix<double> uvw(3, n_baselines_);
    for (size_t i = 0; i < n_baselines_; ++i) {
      uvw(0, i) = 1 + n_time_processed_ + i;
      uvw(1, i) = 2 + n_time_processed_ + i;
      uvw(2, i) = 3 + n_time_processed_ + i;
    }
    DPBuffer buffer;
    buffer.setTime(n_time_processed_ * 5 +
                   3);  // same interval as in updateAveragInfo
    buffer.setData(data);
    buffer.setUVW(uvw);
    casacore::Cube<float> weights(data.shape());
    weights = 1.0;
    buffer.setWeights(weights);
    casacore::Cube<bool> flags(data.shape());
    flags = flag_;
    buffer.setFlags(flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    casacore::Cube<bool> fullResFlags(n_channels_, 1, n_baselines_);
    fullResFlags = flag_;
    buffer.setFullResFlags(fullResFlags);
    getNextStep()->process(buffer);
    ++n_time_processed_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void show(std::ostream&) const override{};
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  size_t n_time_processed_;
  size_t n_time_;
  size_t n_baselines_;
  size_t n_channels_;
  size_t n_correlations_;
  bool flag_;
};

// Class to check result of flagged, unaveraged TestInput run by test1.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int n_time, int n_baselines, int n_channels, int n_correlations)
      : n_time_processed_(0),
        n_time_(n_time),
        n_baselines_(n_baselines),
        n_channels_(n_channels),
        n_correlations_(n_correlations) {}

 private:
  bool process(const DPBuffer& buffer) override {
    casacore::Cube<bool> result(n_correlations_, n_channels_, n_baselines_);
    for (int i = 0; i < n_baselines_; ++i) {
      casacore::Cube<bool> test = buffer.GetCasacoreFlags();

      if (i == 4 || i == 10 || i == 11) {
        for (int j = 0; j < n_channels_; ++j) {
          for (int k = 0; k < n_correlations_; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    BOOST_CHECK(allEQ(buffer.GetCasacoreFlags(), result));
    n_time_processed_++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(static_cast<int>(infoIn.origNChan()), n_channels_);
    BOOST_CHECK_EQUAL(static_cast<int>(infoIn.nchan()), n_channels_);
    BOOST_CHECK_EQUAL(static_cast<int>(infoIn.ntime()), n_time_);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(static_cast<int>(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(static_cast<int>(infoIn.ntimeAvg()), 1);
  }

  int n_time_processed_;
  int n_time_;
  int n_baselines_;
  int n_channels_;
  int n_correlations_;
  bool flag_;
  bool clear_;
};

BOOST_AUTO_TEST_CASE(execute) {
  const int kNTimes = 10;
  const int kNBaselines = 576;
  const int kNChannels = 8;
  const int kNCorrelations = 4;
  const bool kFlag = false;

  auto in = std::make_shared<TestInput>(kNTimes, kNBaselines, kNChannels,
                                        kNCorrelations, kFlag);
  dp3::common::ParameterSet parset;
  parset.add("selection", "4,10,11");
  auto antenna_flagger = std::make_shared<AntennaFlagger>(parset, "");
  auto out = std::make_shared<TestOutput>(kNTimes, kNBaselines, kNChannels,
                                          kNCorrelations);
  dp3::steps::test::Execute({in, antenna_flagger, out});
}

BOOST_AUTO_TEST_SUITE_END()
