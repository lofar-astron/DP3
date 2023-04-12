// tApplyCal.cc: Test program for class ApplyCal
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "../../ApplyCal.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MPosition.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::ApplyCal;
using dp3::steps::Step;

namespace {
const size_t kNTimes = 10;
const double kFirstTime = 4472025740.0;
const int kTimeInterval = 5;
const size_t kNBaselines = 9;
const size_t kNChannels = 32;
const size_t kNCorrelations = 4;
}  // namespace

BOOST_AUTO_TEST_SUITE(applycal)

// Simple class to generate input arrays.
// 9 baselines, 3 antennas, 4 correlations
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput() : process_count_(0) {
    info() = DPInfo(kNCorrelations, kNChannels);
    info().setTimes(kFirstTime, kFirstTime + (kNTimes - 1) * kTimeInterval,
                    kTimeInterval);
    // Fill the baseline stations; use 3 stations.
    // So they are called 00 01 02 10 11 12 20 21 22, etc.

    std::vector<int> ant1(kNBaselines);
    std::vector<int> ant2(kNBaselines);
    int st1 = 0;
    int st2 = 0;
    for (size_t i = 0; i < kNBaselines; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 3) {
        st2 = 0;
        if (++st1 == 3) {
          st1 = 0;
        }
      }
    }
    std::vector<string> antNames{"ant1", "ant2", "ant3", ""};
    // Define their positions (more or less WSRT RT0-3).
    std::vector<casacore::MPosition> antPos(4);
    casacore::Vector<double> vals(3);
    vals[0] = 3828763;
    vals[1] = 442449;
    vals[2] = 5064923;
    antPos[0] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828746;
    vals[1] = 442592;
    vals[2] = 5064924;
    antPos[1] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828729;
    vals[1] = 442735;
    vals[2] = 5064925;
    antPos[2] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828713;
    vals[1] = 442878;
    vals[2] = 5064926;
    antPos[3] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    std::vector<double> antDiam(4, 70.);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(kNChannels, 1000000.);
    std::vector<double> chanFreqs;
    for (size_t i = 0; i < kNChannels; ++i) {
      chanFreqs.push_back(10500000. + i * 1000000.);
    }
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
  }

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (process_count_ == kNTimes) {
      return false;
    }
    buffer->ResizeData(kNBaselines, kNChannels, kNCorrelations);
    buffer->GetData().fill(std::complex<float>(1, 0));
    buffer->ResizeWeights(kNBaselines, kNChannels, kNCorrelations);
    buffer->GetWeights().fill(1.0f);

    buffer->setTime(kFirstTime + process_count_ * kTimeInterval);
    buffer->ResizeFlags(kNBaselines, kNChannels, kNCorrelations);
    buffer->GetFlags().fill(false);
    getNextStep()->process(std::move(buffer));
    ++process_count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

 private:
  int process_count_;
};

// Class to check result of TestInput run by tests.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  enum tests {
    WeightsNotChanged = 1,
    DataNotChanged = 2,
    DataChanged = 4,
    DataEquals = 8,
    WeightEquals = 16
  };
  TestOutput(int enabled_tests)
      : time_step_(0), enabled_tests_(enabled_tests) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Fill data and scale as needed.
    xt::xtensor<std::complex<float>, 3> data(
        {kNBaselines, kNChannels, kNCorrelations}, std::complex<float>(1, 0));

    // The same gain corrections as in tApplyCal_tmp.parmdb
    // Shape: correlation, antenna, time, frequency.
    xt::xtensor<std::complex<float>, 4> gains({4, 3, 2, 2});
    xt::view(gains, 0, xt::all(), xt::all(), xt::all()) = 1.0f;
    xt::view(gains, 1, xt::all(), xt::all(), xt::all()) = 0.0f;
    xt::view(gains, 2, xt::all(), xt::all(), xt::all()) = 0.0f;
    xt::view(gains, 3, xt::all(), xt::all(), xt::all()) = 1.0f;
    // ant2
    gains(0, 1, 0, 0) = 2.0f;
    gains(3, 1, 0, 0) = 3.0f;
    gains(0, 1, 1, 1) = std::complex<float>(3.0f, 4.0f);
    // ant3
    gains(2, 2, 0, 1) = 0.5f;

    if (enabled_tests_ & DataEquals) {
      const int time = time_step_ / (kNTimes / 2);

      for (size_t baseline = 0; baseline < kNBaselines; ++baseline) {
        const int antenna1 = info().getAnt1()[baseline];
        const int antenna2 = info().getAnt2()[baseline];

        for (size_t channel = 0; channel < kNChannels; ++channel) {
          const size_t gain_channel = channel / (kNChannels / 2);

          for (size_t corr = 0; corr < kNCorrelations; ++corr) {
            data(baseline, channel, corr) /=
                gains(corr / 2 * 3, antenna1, time, gain_channel) *
                std::conj(gains(corr % 2 * 3, antenna2, time, gain_channel));
          }

          if (antenna2 == 2 && time == 0 && gain_channel == 1) {
            data(baseline, channel, 0) -= 0.5f * data(baseline, channel, 1);
            data(baseline, channel, 2) -= 0.5f * data(baseline, channel, 3);
          }

          if (antenna1 == 2 && time == 0 && gain_channel == 1) {
            data(baseline, channel, 0) -= 0.5f * data(baseline, channel, 2);
            data(baseline, channel, 1) -= 0.5f * data(baseline, channel, 3);
          }
        }
      }
    }

    if (enabled_tests_ & WeightEquals) {
      BOOST_CHECK_CLOSE(buffer->GetWeights()(1, 0, 0), 4.0, 1.0e-3);
      BOOST_CHECK_CLOSE(buffer->GetWeights()(1, 0, 1), 9.0, 1.0e-3);
      BOOST_CHECK_CLOSE(buffer->GetWeights()(1, 0, 2), 4.0, 1.0e-3);
      BOOST_CHECK_CLOSE(buffer->GetWeights()(1, 0, 3), 9.0, 1.0e-3);
      BOOST_CHECK_CLOSE(buffer->GetWeights()(5, 31, 0), 0.8, 1.0e-3);
    }

    if (enabled_tests_ & (DataEquals | DataNotChanged)) {
      BOOST_CHECK(xt::allclose(xt::real(buffer->GetData()), xt::real(data),
                               1.0e-7, 1.0e-10));
      BOOST_CHECK(xt::allclose(xt::imag(buffer->GetData()), xt::imag(data),
                               1.0e-7, 1.0e-10));
    }
    if (enabled_tests_ & DataChanged) {
      BOOST_CHECK(!xt::allclose(xt::real(buffer->GetData()), xt::real(data),
                                1.0e-7, 1.0e-10));
    }
    if (enabled_tests_ & WeightsNotChanged) {
      BOOST_CHECK(xt::allclose(buffer->GetWeights(),
                               xt::ones<float>(data.shape()), 1.0e-6, 1.0e-10));
    }
    time_step_++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(kNChannels, infoIn.origNChan());
    BOOST_CHECK_EQUAL(kNChannels, infoIn.nchan());
    BOOST_CHECK_EQUAL(kNTimes, infoIn.ntime());
    BOOST_CHECK_EQUAL(kTimeInterval, infoIn.timeInterval());
    BOOST_CHECK_EQUAL(kNBaselines, infoIn.nbaselines());
  }

 private:
  int time_step_;
  int enabled_tests_;
};

// Test clock + tec, and test two ApplyCals in sequence
BOOST_AUTO_TEST_CASE(test_clock_and_tec) {
  auto input = std::make_shared<TestInput>();

  dp3::common::ParameterSet parset1;
  parset1.add("correction", "tec");
  parset1.add("parmdb", "tApplyCal_tmp.parmdb");
  parset1.add("timeslotsperparmupdate", "5");
  parset1.add("updateweights", "true");
  auto apply_cal1 = std::make_shared<ApplyCal>(parset1, "");

  dp3::common::ParameterSet parset2;
  parset2.add("correction", "clock");
  parset2.add("parmdb", "tApplyCal_tmp.parmdb");
  parset2.add("timeslotsperparmupdate", "5");
  parset2.add("updateweights", "true");
  auto apply_cal2 = std::make_shared<ApplyCal>(parset2, "");

  dp3::common::ParameterSet parset3;
  parset3.add("correction", "commonscalarphase");
  parset3.add("parmdb", "tApplyCal_tmp.parmdb");
  parset3.add("timeslotsperparmupdate", "1");
  parset3.add("udpateweights", "true");
  auto apply_cal3 = std::make_shared<ApplyCal>(parset3, "");

  auto output = std::make_shared<TestOutput>(TestOutput::DataChanged |
                                             TestOutput::WeightsNotChanged);

  dp3::steps::test::Execute(
      {input, apply_cal1, apply_cal2, apply_cal3, output});
}

BOOST_AUTO_TEST_CASE(test_gain) {
  auto input = std::make_shared<TestInput>();

  dp3::common::ParameterSet parset1;
  parset1.add("correction", "gain");
  parset1.add("parmdb", "tApplyCal_tmp.parmdb");
  parset1.add("timeslotsperparmupdate", "5");
  parset1.add("updateweights", "true");
  auto apply_cal = std::make_shared<ApplyCal>(parset1, "");

  auto output = std::make_shared<TestOutput>(TestOutput::DataEquals);

  dp3::steps::test::Execute({input, apply_cal, output});
}

BOOST_AUTO_TEST_SUITE_END()
