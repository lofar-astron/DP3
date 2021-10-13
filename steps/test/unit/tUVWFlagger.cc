// tUVWFlagger.cc: Test program for class UVWFlagger
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "tStepCommon.h"

#include "../../UVWFlagger.h"
#include "../../InputStep.h"

#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"

#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::Step;
using std::vector;

namespace {

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::InputStep {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {
    info().init(ncorr, 0, nchan, ntime, 0., 5., string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(nbl);
    vector<int> ant2(nbl);
    int st1 = 0;
    int st2 = 0;
    for (int i = 0; i < nbl; ++i) {
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
    std::vector<double> chanWidth(nchan, 1000000.);
    std::vector<double> chanFreqs;
    for (int i = 0; i < nchan; i++) {
      chanFreqs.push_back(10500000. + i * 1000000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }

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
    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime(itsCount * 30 + 4472025740.0);
    buf.setData(data);
    buf.setUVW(uvw);
    casacore::Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags(flags);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of flagged, unaveraged TestInput run by test1.
class TestOutput : public Step {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Flag where u,v,w matches intervals given in test1.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    result = false;
    for (int i = 0; i < itsNBl; ++i) {
      double u = 1 + i + itsCount;
      double v = 2 + i + itsCount;
      double w = 3 + i + itsCount;
      double uv = sqrt(u * u + v * v);
      if ((uv > 5.5 && uv < 8.5) || (u > 20.5 && u < 23.5) ||
          (u > 31.5 && u < 40.5) || (v > 11.5 && v < 14.5) || w < 3.5 ||
          w > 44.5) {
        for (int j = 0; j < itsNChan; ++j) {
          for (int k = 0; k < itsNCorr; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    BOOST_CHECK(allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of flagged, unaveraged TestInput run by test2.
class TestOutput2 : public Step {
 public:
  TestOutput2(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Flag where u,v,w matches intervals given in test1.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    result = false;
    for (int i = 0; i < itsNBl; ++i) {
      for (int j = 0; j < itsNChan; ++j) {
        double wavel = 2.99792458e+08 / (10.5e6 + j * 1e6);
        double u = (1 + i + itsCount) / wavel;
        double v = (2 + i + itsCount) / wavel;
        double w = (3 + i + itsCount) / wavel;
        double uv = sqrt(u * u + v * v);
        if ((uv > 0.2 && uv < 0.31) || (u > 1.55 && u < 1.485) ||
            (u > 0.752 && u < 0.862) || (v > 0.42 && v < 0.53) || w < 0.12 ||
            w > 1.63) {
          for (int k = 0; k < itsNCorr; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    BOOST_CHECK(allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of flagged, unaveraged TestInput run by test3.
class TestOutput3 : public Step {
 public:
  TestOutput3(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {
    BOOST_CHECK_EQUAL(ntime, 2);
    BOOST_CHECK_EQUAL(nbl, 16);
  }

 private:
  virtual bool process(const DPBuffer& buf) {
    // These are the UVW coordinates as calculated by UVWFlagger for the
    // station positions and times defined in TestInput and phase center
    // defined in test3.
    double uvwvals[] = {
        0,         0,        0,        0.423756,  -127.372, 67.1947,
        0.847513,  -254.744, 134.389,  0.277918,  -382.015, 201.531,
        -0.423756, 127.372,  -67.1947, 0,         0,        0,
        0.423756,  -127.372, 67.1947,  -0.145838, -254.642, 134.336,
        -0.847513, 254.744,  -134.389, -0.423756, 127.372,  -67.1947,
        0,         0,        0,        -0.569594, -127.27,  67.1417,
        -0.277918, 382.015,  -201.531, 0.145838,  254.642,  -134.336,
        0.569594,  127.27,   -67.1417, 0,         0,        0,
        0,         0,        0,        0.738788,  -127.371, 67.1942,
        1.47758,   -254.742, 134.388,  1.22276,   -382.013, 201.53,
        -0.738788, 127.371,  -67.1942, 0,         0,        0,
        0.738788,  -127.371, 67.1942,  0.483976,  -254.642, 134.336,
        -1.47758,  254.742,  -134.388, -0.738788, 127.371,  -67.1942,
        0,         0,        0,        -0.254812, -127.271, 67.1421,
        -1.22276,  382.013,  -201.53,  -0.483976, 254.642,  -134.336,
        0.254812,  127.271,  -67.1421, 0,         0,        0};
    casacore::Cube<double> uvws(casacore::IPosition(3, 3, 16, 2), uvwvals,
                                casacore::SHARE);
    // Flag where u,v,w matches intervals given in test3.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    result = false;
    for (int i = 0; i < itsNBl; ++i) {
      double u = uvws(0, i, itsCount);
      double v = uvws(1, i, itsCount);
      double w = uvws(2, i, itsCount);
      double uv = sqrt(u * u + v * v);
      if ((uv > 5.5 && uv < 8.5) || (u > 20.5 && u < 23.5) ||
          (u > 31.5 && u < 40.5) || (v > 11.5 && v < 14.5) || w < 3.5 ||
          w > 44.5) {
        for (int j = 0; j < itsNChan; ++j) {
          for (int k = 0; k < itsNCorr; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    BOOST_CHECK(allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Test flagging a few baselines on UV in m.
void test1(int ntime, int nbl, int nchan, int ncorr) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  Step::ShPtr step1(in);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add("vmrange", "[11.5..14.5]");
  parset.add("wmmax", "44.5");
  parset.add("wmmin", "3.5");
  Step::ShPtr step2(new dp3::steps::UVWFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test flagging a few baselines on UV in wavelengths.
void test2(int ntime, int nbl, int nchan, int ncorr) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  Step::ShPtr step1(in);
  dp3::common::ParameterSet parset;
  parset.add("uvlambdarange", "[0.2..0.31]");
  parset.add("ulambdarange", "[1.55..1.485, 0.807+-0.055]");
  parset.add("vlambdarange", "[0.42..0.53]");
  parset.add("wlambdamax", "1.63");
  parset.add("wlambdamin", "0.12");
  Step::ShPtr step2(new dp3::steps::UVWFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test flagging a few baselines on UV in m with a different phase center.
void test3(int ntime, int nbl, int nchan, int ncorr) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  Step::ShPtr step1(in);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add("vmrange", "[11.5..14.5]");
  parset.add("wmmax", "44.5");
  parset.add("wmmin", "3.5");
  parset.add("phasecenter", "[-1.92653768rad, 1.09220917rad, j2000]");
  Step::ShPtr step2(new dp3::steps::UVWFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput3(ntime, nbl, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

}  // namespace

BOOST_AUTO_TEST_SUITE(uvwflagger)

BOOST_AUTO_TEST_CASE(testuvwflagger1) { test1(10, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(testuvwflagger2) { test1(100, 105, 32, 4); }

BOOST_AUTO_TEST_CASE(testuvwflagger3) { test2(2, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(testuvwflagger4) { test2(2, 36, 16, 2); }

BOOST_AUTO_TEST_CASE(testuvwflagger5) { test2(10, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(testuvwflagger6) { test2(100, 105, 32, 4); }

BOOST_AUTO_TEST_CASE(testuvwflagger7) { test3(2, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(sun_as_phase_center) {
  auto in = std::make_shared<TestInput>(1, 1, 1, 1);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("phasecenter", "Sun");
  BOOST_CHECK_NO_THROW(
      std::make_shared<dp3::steps::UVWFlagger>(in.get(), parset, ""));
}

BOOST_AUTO_TEST_SUITE_END()
