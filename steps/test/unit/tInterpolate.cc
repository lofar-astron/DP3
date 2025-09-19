// tInterpolate.cc: Test program for class Interpolate
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Interpolate.h"

#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPInfo;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(interpolate)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(std::size_t num_times, std::size_t num_antennas,
            std::size_t num_channels, std::size_t num_correlations, bool flag)
      : process_count_(0),
        num_times_(num_times),
        num_baselines_(num_antennas * (num_antennas + 1) / 2),
        num_channels_(num_channels),
        num_correlations_(num_correlations),
        flag_(flag) {
    using casacore::MPosition;
    using casacore::Quantum;

    GetWritableInfoOut() = DPInfo(num_correlations_, num_channels_);
    GetWritableInfoOut().setTimes(100.0, 100.0 + (num_times_ - 1) * 5.0, 5.0);

    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> antenna_1(num_baselines_);
    std::vector<int> antenna_2(num_baselines_);
    int st1 = 0;
    int st2 = 0;
    for (std::size_t i = 0; i < num_baselines_; ++i) {
      antenna_1[i] = st1;
      antenna_2[i] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antenna_names{"rs01.s01", "rs02.s01", "cs01.s01",
                                           "cs01.s02"};
    // Define their positions (more or less WSRT RT0-3).
    const casacore::Vector<double> vals0(
        std::vector<int>{3828763, 442449, 5064923});
    const casacore::Vector<double> vals1(
        std::vector<int>{3828746, 442592, 5064924});
    const casacore::Vector<double> vals2(
        std::vector<int>{3828729, 442735, 5064925});
    const casacore::Vector<double> vals3(
        std::vector<int>{3828713, 442878, 5064926});
    const std::vector<MPosition> antenna_positions{
        MPosition(Quantum<casacore::Vector<double>>(vals0, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals1, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals2, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals3, "m"),
                  MPosition::ITRF)};
    const std::vector<double> antenna_diameters(4.0, 70.0);
    GetWritableInfoOut().setAntennas(antenna_names, antenna_diameters,
                                     antenna_positions, antenna_1, antenna_2);
    // Define the frequencies.
    std::vector<double> channel_frequencies;
    std::vector<double> channel_widths(num_channels, 100000.0);
    for (std::size_t i = 0; i < num_channels; ++i) {
      channel_frequencies.push_back(1050000.0 + i * 100000.0);
    }
    GetWritableInfoOut().setChannels(std::move(channel_frequencies),
                                     std::move(channel_widths));
  }

 private:
  bool process(std::unique_ptr<dp3::base::DPBuffer>) override {
    // Stop when all times are done.
    if (process_count_ == num_times_) {
      return false;
    }
    const std::array<std::size_t, 3> data_shape{num_baselines_, num_channels_,
                                                num_correlations_};
    xt::xtensor<std::complex<float>, 3> data(data_shape);
    for (std::size_t i = 0; i < data.size(); ++i) {
      data.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (process_count_ == 5) {
      data += std::complex<float>(10.0, 10.0);
    }
    auto buffer = std::make_unique<dp3::base::DPBuffer>();
    buffer->SetTime(process_count_ * 5.0 +
                    2.0);  // same interval as in updateAveragInfo
    buffer->GetData().resize(data_shape);
    buffer->GetData() = data;
    xt::xtensor<float, 3> weights(data_shape, 1.0);
    buffer->GetWeights().resize(data_shape);
    buffer->GetWeights() = weights;
    xt::xtensor<bool, 3> flags(data_shape, flag_);
    buffer->GetFlags().resize(data_shape);
    buffer->GetFlags() = flags;
    getNextStep()->process(std::move(buffer));
    ++process_count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  std::size_t process_count_, num_times_, num_baselines_, num_channels_,
      num_correlations_;
  bool flag_;
};

// Class to check result.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(std::size_t num_times, std::size_t num_antennas,
             std::size_t num_channels, std::size_t num_correlations)
      : process_count_(0),
        num_times_(num_times),
        num_baselines_(num_antennas * (num_antennas + 1) / 2),
        num_channels_(num_channels),
        num_correlations_(num_correlations) {}

 private:
  bool process(
      const std::unique_ptr<dp3::base::DPBuffer> input_buffer) override {
    // Fill expected result in similar way as TestInput.
    const std::array<std::size_t, 3> data_shape{num_baselines_, num_channels_,
                                                num_correlations_};
    xt::xtensor<std::complex<float>, 3> expected_result(data_shape);
    for (std::size_t i = 0; i < expected_result.size(); ++i) {
      expected_result.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (process_count_ == 5) {
      expected_result += std::complex<float>(10.0, 10.0);
    }
    // Check the result.
    BOOST_CHECK(xt::allclose(input_buffer->GetData(), expected_result, 1e-10));
    BOOST_CHECK(
        xt::allclose(input_buffer->GetTime(), 2 + 5.0 * process_count_));
    ++process_count_;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(info.origNChan(), num_channels_);
    BOOST_CHECK_EQUAL(info.nchan(), num_channels_);
    BOOST_CHECK_EQUAL(info.ntime(), num_times_);
    BOOST_CHECK_EQUAL(info.firstTime(), 100.0);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5.0);
    BOOST_CHECK_EQUAL(info.nchanAvg(), 1);
    BOOST_CHECK_EQUAL(info.ntimeAvg(), 1);
    BOOST_CHECK_EQUAL(info.chanFreqs().size(), num_channels_);
    BOOST_CHECK_EQUAL(info.chanWidths().size(), num_channels_);
    BOOST_CHECK(info.msName().empty());
  }

  std::size_t process_count_, num_times_, num_baselines_, num_channels_,
      num_correlations_;
};

BOOST_AUTO_TEST_CASE(test1) {
  const int kNTime = 10;
  const int kNAnt = 2;
  const int kNChan = 32;
  const int kNCorr = 4;
  const bool kFlag = false;

  auto step_1 =
      std::make_shared<TestInput>(kNTime, kNAnt, kNChan, kNCorr, kFlag);
  dp3::common::ParameterSet parset;
  parset.add("windowsize", "9");
  auto step_2 = std::make_shared<dp3::steps::Interpolate>(parset, "");
  auto step_3 = std::make_shared<TestOutput>(kNTime, kNAnt, kNChan, kNCorr);
  dp3::steps::test::Execute({step_1, step_2, step_3});
}

BOOST_AUTO_TEST_SUITE_END()
