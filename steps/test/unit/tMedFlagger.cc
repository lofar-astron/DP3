// tMedFlagger.cc: Test program for class MedFlagger
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "tStepCommon.h"
#include "../../MedFlagger.h"
#include "../../InputStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::MedFlagger;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(medflagger)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public InputStep {
 public:
  TestInput(int ntime, int nant, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 10 + itsCount * 6);
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
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    // Use timeInterval=5
    info().init(itsNCorr, 0, itsNChan, itsNTime, 100, 5, string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(itsNBl);
    vector<int> ant2(itsNBl);
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
    vector<string> antNames{"rs01.s01", "rs02.s01", "cs01.s01", "cs01.s02"};
    // Define their positions (more or less WSRT RT0-3).
    vector<casacore::MPosition> antPos(4);
    vector<double> vals(3);
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
    vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNChan, 100000);
    for (int i = 0; i < itsNChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result.
class TestOutput : public Step {
 public:
  TestOutput(int ntime, int nant, int nchan, int ncorr, bool flag,
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
  virtual bool process(const DPBuffer& buf) {
    // Fill expected result in similar way as TestInput.
    casacore::Cube<casacore::Complex> result(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(result.size()); ++i) {
      result.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 10 + itsCount * 6);
    }
    // Check the result.
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-10));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-10));
    // Check the flags.
    // If autocorrs are used, only the last channel is flagged, but the first
    // channel also for the first time stamp. Thus is only true for a limited
    // nr of baselines (thus do not use nant>2 in test2 with flag=false).
    // If short baselines are used, bl 2,3,7,8,12,13 are not flagged.
    // The others have length 0 or 144.
    casacore::Cube<bool> expFlag(itsNCorr, itsNChan, itsNBl);
    expFlag = itsFlag;
    if (itsUseAutoCorr) {
      for (int i = 0; i < itsNBl; ++i) {
        if (!itsShortBL ||
            !(i == 2 || i == 3 || i == 7 || i == 8 || i == 12 || i == 13)) {
          for (int j = 0; j < itsNCorr; ++j) {
            expFlag(j, 0, i) = itsFlag || itsCount == 0;
            expFlag(j, itsNChan - 1, i) = true;
          }
        }
      }
    }
    BOOST_CHECK(allEQ(buf.getFlags(), expFlag));
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
    BOOST_CHECK_EQUAL(info.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(info.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(info.ntimeAvg()), 1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag, itsUseAutoCorr, itsShortBL;
};

// Test simple flagging with or without preflagged points.
void test1(int ntime, int nant, int nchan, int ncorr, bool flag, int threshold,
           bool shortbl) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("freqwindow", "1");
  parset.add("timewindow", "1");
  parset.add("threshold", std::to_string(threshold));
  if (shortbl) {
    parset.add("blmin", "0");
    parset.add("blmax", "145");
  }
  Step::ShPtr step2(new MedFlagger(in, parset, ""));
  Step::ShPtr step3(
      new TestOutput(ntime, nant, nchan, ncorr, flag, false, shortbl));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test applyautocorr flagging with or without preflagged points.
void test2(int ntime, int nant, int nchan, int ncorr, bool flag, int threshold,
           bool shortbl) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("freqwindow", "3");
  parset.add("timewindow", "min(1,max(1,bl))");
  parset.add("threshold", std::to_string(threshold));
  parset.add("applyautocorr", "True");
  if (shortbl) {
    parset.add("blmax", "145");
  }
  Step::ShPtr step2(new MedFlagger(in, parset, ""));
  Step::ShPtr step3(
      new TestOutput(ntime, nant, nchan, ncorr, flag, true, shortbl));
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_DATA_TEST_CASE(test_medflagger_1,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(10, 2, 32, 4, false, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_medflagger_2,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(10, 5, 32, 4, true, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_medflagger_3,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test1(4, 2, 8, 4, false, 100, shortbl);
}

BOOST_DATA_TEST_CASE(test_medflagger_4,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test2(10, 5, 32, 4, true, 1, shortbl);
}

BOOST_DATA_TEST_CASE(test_medflagger_5,
                     boost::unit_test::data::make({true, false}), shortbl) {
  test2(4, 2, 8, 4, false, 100, shortbl);
}

BOOST_AUTO_TEST_SUITE_END()
