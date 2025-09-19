// tAOFlaggerStep.cc: Test program for class AOFlaggerStep
//
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <boost/test/unit_test.hpp>
#include <xtensor/xcomplex.hpp>

#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../../AOFlaggerStep.h"
#include "../../../base/test/LoggerFixture.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::AOFlaggerStep;
using dp3::steps::Step;

using dp3::common::ParameterSet;

using casacore::MPosition;
using casacore::Quantum;

namespace {

static const float kFirstTime = 100.0;
static const float kTimeStep = 5.0;
static const size_t kOutlierIndex = 5;

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(int ntime, int nant, int nchan, int ncorr, bool flag)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nant * (nant + 1) / 2),
        n_channels_(nchan),
        n_correlations_(ncorr),
        flag_(flag) {
    GetWritableInfoOut() = DPInfo(n_correlations_, n_channels_);
    GetWritableInfoOut().setTimes(
        kFirstTime, kFirstTime + (n_times_ - 1) * kTimeStep, kTimeStep);

    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(n_baselines_);
    std::vector<int> ant2(n_baselines_);
    int st1 = 0;
    int st2 = 0;
    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      ant1[bl] = st1;
      ant2[bl] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antNames(4);
    antNames[0] = "rs01.s01";
    antNames[1] = "rs02.s01";
    antNames[2] = "cs01.s01";
    antNames[3] = "cs01.s02";
    // Define their positions (more or less WSRT RT0-3).
    std::vector<casacore::MPosition> antPos(4);
    casacore::Vector<double> vals(3);
    vals[0] = 3828763;
    vals[1] = 442449;
    vals[2] = 5064923;
    antPos[0] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828746;
    vals[1] = 442592;
    vals[2] = 5064924;
    antPos[1] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828729;
    vals[1] = 442735;
    vals[2] = 5064925;
    antPos[2] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828713;
    vals[1] = 442878;
    vals[2] = 5064926;
    antPos[3] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    std::vector<double> antDiam(4, 70.);
    GetWritableInfoOut().setAntennas(antNames, antDiam, antPos, ant1, ant2);

    // Define the frequencies.
    std::vector<double> chanFreqs(nchan);
    std::vector<double> chanWidth(nchan, 100000.);
    std::iota(chanFreqs.begin(), chanFreqs.end(), 1050000);
    GetWritableInfoOut().setChannels(std::move(chanFreqs),
                                     std::move(chanWidth));
  }

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (count_ == n_times_) {
      return false;
    }
    const std::array<std::size_t, 3> shape{n_baselines_, n_channels_,
                                           n_correlations_};
    buffer->GetData().resize(shape);

    buffer->GetData().fill(std::complex<float>(1.6, 0.9));
    if (count_ == kOutlierIndex) {
      buffer->GetData() += std::complex<float>(10., 10.);
    }
    buffer->SetTime(count_ * kTimeStep +
                    kFirstTime);  // same interval as in updateAveragInfo
    buffer->GetWeights().resize(shape);
    buffer->GetWeights().fill(1.);
    buffer->GetFlags().resize(shape);
    buffer->GetFlags().fill(flag_);
    getNextStep()->process(std::move(buffer));
    ++count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  size_t count_, n_times_, n_baselines_, n_channels_, n_correlations_;
  bool flag_;
};

// Class to check result.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, int nant, int nchan, int ncorr)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nant * (nant + 1) / 2),
        n_channels_(nchan),
        n_correlations_(ncorr) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Fill expected result in similar way as TestInput.
    xt::xtensor<std::complex<float>, 3> result(
        {n_baselines_, n_channels_, n_correlations_});
    for (int i = 0; i < int(result.size()); ++i) {
      result.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (count_ == kOutlierIndex) {
      result += std::complex<float>(10.0, 10.0);
    }
    // Check the result.
    BOOST_CHECK(xt::allclose(buffer->GetData(), result, 1.0e-10));
    BOOST_CHECK_CLOSE(buffer->GetTime(), kFirstTime + kTimeStep * count_,
                      1.0e-8);
    ++count_;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK(info.origNChan() == n_channels_);
    BOOST_CHECK(info.nchan() == n_channels_);
    BOOST_CHECK(info.ntime() == n_times_);
    BOOST_CHECK(info.firstTime() == kFirstTime);
    BOOST_CHECK(info.timeInterval() == kTimeStep);
    BOOST_CHECK(int(info.nchanAvg()) == 1);
    BOOST_CHECK(int(info.ntimeAvg()) == 1);
    BOOST_CHECK(info.chanFreqs().size() == n_channels_);
    BOOST_CHECK(info.chanWidths().size() == n_channels_);
    BOOST_CHECK(info.msName().empty());
  }

  size_t count_;
  size_t n_times_, n_baselines_, n_channels_, n_correlations_;
};

// Test simple flagging with or without preflagged points.
void test1(int ntime, int nant, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("timewindow", "1");
  Step::ShPtr step2(new AOFlaggerStep(parset, ""));
  Step::ShPtr step3(new TestOutput(ntime, nant, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test applyautocorr flagging with or without preflagged points.
void test2(int ntime, int nant, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("timewindow", "4");
  parset.add("overlapmax", "1");
  Step::ShPtr step2(new AOFlaggerStep(parset, ""));
  Step::ShPtr step3(new TestOutput(ntime, nant, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

}  // namespace

BOOST_AUTO_TEST_SUITE(
    aoflaggerstep, *boost::unit_test::fixture<dp3::base::test::LoggerFixture>())

BOOST_AUTO_TEST_CASE(legacy_test1) {
  for (unsigned int i = 0; i < 2; ++i) {
    test1(10, 2, 32, 4, false);
    test1(10, 5, 32, 4, true);
  }
}

BOOST_AUTO_TEST_CASE(legacy_test2) {
  for (unsigned int i = 0; i < 2; ++i) {
    test2(4, 2, 8, 4, false);
    test2(10, 5, 32, 4, true);
    test2(8, 2, 8, 4, false);
    test2(14, 2, 8, 4, false);
  }
}

BOOST_AUTO_TEST_SUITE_END()
