// tMadFlagger.cc: Test program for class MadFlagger
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../MadFlagger.h"

#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MPosition.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::MadFlagger;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(madflagger)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(size_t ntime, size_t nant, size_t nchan, size_t ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }

    buffer->setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo

    const std::array<size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    buffer->ResizeData(shape);
    for (size_t i = 0; i < buffer->GetData().size(); ++i) {
      buffer->GetData().data()[i] = std::complex<float>(
          i + itsCount * 10.0f, i - 10.0f + itsCount * 6.0f);
    }

    buffer->ResizeWeights(shape);
    buffer->GetWeights().fill(1.0);

    buffer->ResizeFlags(shape);
    buffer->GetFlags().fill(itsFlag);

    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    info() = DPInfo(itsNCorr, itsNChan);
    // Use timeInterval=5
    info().setTimes(100.0, 100.0 + (itsNTime - 1) * 5.0, 5.0);
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(itsNBl);
    std::vector<int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (size_t bl = 0; bl < itsNBl; ++bl) {
      ant1[bl] = st1;
      ant2[bl] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antNames{"rs01.s01", "rs02.s01", "cs01.s01",
                                      "cs01.s02"};
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
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNChan, 100000);
    for (size_t chan = 0; chan < itsNChan; chan++) {
      chanFreqs.push_back(1050000. + chan * 100000.);
    }
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
  }

  size_t itsCount;
  size_t itsNTime;
  size_t itsNBl;
  size_t itsNChan;
  size_t itsNCorr;
  bool itsFlag;
};

// Class to check result.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(size_t ntime, size_t nant, size_t nchan, size_t ncorr, bool flag,
             bool useAutoCorr, bool shortbl)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag),
        itsUseAutoCorr(useAutoCorr),
        itsShortBL(shortbl) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    std::array<size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    // Fill expected result in similar way as TestInput.
    xt::xtensor<std::complex<float>, 3> result_data(shape, 0.0f);
    for (size_t i = 0; i < result_data.size(); ++i) {
      result_data.data()[i] = std::complex<float>(i + itsCount * 10.0f,
                                                  i - 10.0f + itsCount * 6.0f);
    }
    // Check the result.
    BOOST_CHECK(xt::allclose(buffer->GetData(), result_data));

    // Check the flags.
    // If autocorrs are used, only the last channel is flagged, but the first
    // channel also for the first time stamp. Thus is only true for a limited
    // nr of baselines (thus do not use nant>2 in test2 with flag=false).
    // If short baselines are used, bl 2,3,7,8,12,13 are not flagged.
    // The others have length 0 or 144.
    xt::xtensor<bool, 3> result_flags(shape, itsFlag);
    if (itsUseAutoCorr) {
      for (size_t bl = 0; bl < itsNBl; ++bl) {
        if (!itsShortBL || !(bl == 2 || bl == 3 || bl == 7 || bl == 8 ||
                             bl == 12 || bl == 13)) {
          for (size_t corr = 0; corr < itsNCorr; ++corr) {
            result_flags(bl, 0, corr) = itsFlag | (itsCount == 0);
            result_flags(bl, itsNChan - 1, corr) = true;
          }
        }
      }
    }
    BOOST_CHECK_EQUAL(buffer->GetFlags(), result_flags);
    BOOST_CHECK_CLOSE(buffer->getTime(), 2 + 5. * itsCount, 1.0e-3);

    ++itsCount;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(size_t(info.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(size_t(info.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(size_t(info.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5);
    BOOST_CHECK_EQUAL(size_t(info.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(size_t(info.ntimeAvg()), 1);
  }

  size_t itsCount;
  size_t itsNTime;
  size_t itsNBl;
  size_t itsNChan;
  size_t itsNCorr;
  bool itsFlag;
  bool itsUseAutoCorr;
  bool itsShortBL;
};

// Test simple flagging with or without preflagged points.
void test1(size_t ntime, size_t nant, size_t nchan, size_t ncorr, bool flag,
           size_t threshold, bool shortbl) {
  // Create the steps.
  auto in = std::make_shared<TestInput>(ntime, nant, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqwindow", "1");
  parset.add("timewindow", "1");
  parset.add("threshold", std::to_string(threshold));
  if (shortbl) {
    parset.add("blmin", "0");
    parset.add("blmax", "145");
  }
  auto mad_flagger = std::make_shared<MadFlagger>(parset, "");
  auto out = std::make_shared<TestOutput>(ntime, nant, nchan, ncorr, flag,
                                          false, shortbl);
  dp3::steps::test::Execute({in, mad_flagger, out});
}

// Test applyautocorr flagging with or without preflagged points.
void test2(size_t ntime, size_t nant, size_t nchan, size_t ncorr, bool flag,
           size_t threshold, bool shortbl) {
  // Create the steps.
  auto in = std::make_shared<TestInput>(ntime, nant, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqwindow", "3");
  parset.add("timewindow", "min(1,max(1,bl))");
  parset.add("threshold", std::to_string(threshold));
  parset.add("applyautocorr", "True");
  if (shortbl) {
    parset.add("blmax", "145");
  }
  auto mad_flagger = std::make_shared<MadFlagger>(parset, "");
  auto out = std::make_shared<TestOutput>(ntime, nant, nchan, ncorr, flag, true,
                                          shortbl);
  dp3::steps::test::Execute({in, mad_flagger, out});
}

BOOST_DATA_TEST_CASE(test_madflagger_1,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(10, 2, 32, 4, false, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_madflagger_2,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(10, 5, 32, 4, true, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_madflagger_3,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(4, 2, 8, 4, false, 100, shortbl);
}

BOOST_DATA_TEST_CASE(test_madflagger_4,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test2(10, 5, 32, 4, true, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_madflagger_5,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test2(4, 2, 8, 4, false, 100, shortbl);
}

BOOST_AUTO_TEST_SUITE_END()
