// tInterpolate.cc: Test program for class Interpolate
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Interpolate.h"

#include "tStepCommon.h"
#include "../../InputStep.h"
#include "../../../base/DP3.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(interpolate)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::InputStep {
 public:
  TestInput(int ntime, int nant, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {
    using casacore::MPosition;
    using casacore::Quantum;

    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(itsNBl);
    std::vector<int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (int i = 0; i < itsNBl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
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
    const std::vector<double> vals0{3828763, 442449, 5064923};
    const std::vector<double> vals1{3828746, 442592, 5064924};
    const std::vector<double> vals2{3828729, 442735, 5064925};
    const std::vector<double> vals3{3828713, 442878, 5064926};
    const std::vector<MPosition> antPos{
        MPosition(Quantum<casacore::Vector<double>>(vals0, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals1, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals2, "m"),
                  MPosition::ITRF),
        MPosition(Quantum<casacore::Vector<double>>(vals3, "m"),
                  MPosition::ITRF)};
    const std::vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(nchan, 100000.);
    for (int i = 0; i < nchan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    casacore::Cube<std::complex<float>> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (itsCount == 5) {
      data += std::complex<float>(10., 10.);
    }
    DPBuffer buf;
    buf.setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags(flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    casacore::Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags(fullResFlags);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {
    // Use startchan=0 and timeInterval=5
    info().init(itsNCorr, 0, itsNChan, itsNTime, 100, 5, string(), string());
  }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result.
class TestOutput : public Step {
 public:
  TestOutput(int ntime, int nant, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Fill expected result in similar way as TestInput.
    casacore::Cube<casacore::Complex> result(itsNCorr, itsNChan, itsNBl);
    for (std::size_t i = 0; i < result.size(); ++i) {
      result.data()[i] = casacore::Complex(1.6, 0.9);
    }
    if (itsCount == 5) {
      result += casacore::Complex(10., 10.);
    }
    // Check the result.
    BOOST_CHECK(casacore::allNear(buf.getData(), result, 1e-10));
    BOOST_CHECK(casacore::near(buf.getTime(), 2 + 5. * itsCount));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(int(info.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(info.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(info.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(info.startTime(), 100);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(info.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(info.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(int(info.chanFreqs().size()), itsNChan);
    BOOST_CHECK_EQUAL(int(info.chanWidths().size()), itsNChan);
    BOOST_CHECK(info.msName().empty());
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

BOOST_AUTO_TEST_CASE(test1) {
  const int kNTime = 10;
  const int kNAnt = 2;
  const int kNChan = 32;
  const int kNCorr = 4;
  const bool kFlag = false;

  auto step1 =
      std::make_shared<TestInput>(kNTime, kNAnt, kNChan, kNCorr, kFlag);
  dp3::common::ParameterSet parset;
  parset.add("windowsize", "9");
  auto step2 =
      std::make_shared<dp3::steps::Interpolate>(step1.get(), parset, "");
  auto step3 = std::make_shared<TestOutput>(kNTime, kNAnt, kNChan, kNCorr);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_SUITE_END()
