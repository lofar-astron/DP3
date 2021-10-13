// tPreFlagger.cc: Test program for class PreFlagger
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../PreFlagger.h"
#include "../../Counter.h"
#include "../../InputStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Counter;
using dp3::steps::InputStep;
using dp3::steps::PreFlagger;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(preflagger)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public InputStep {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {
    // Define start time 0.5 (= 3 - 0.5*5) and time interval 5.
    info().init(ncorr, 0, nchan, ntime, 0.5, 5., string(), string());
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
    casacore::Vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 100000.);
    std::vector<double> chanFreqs;
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
    buf.setTime(itsCount * 5 + 3);  // same interval as in updateAveragInfo
    buf.setData(data);
    buf.setUVW(uvw);
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
  virtual void updateInfo(const DPInfo&) {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of flagged, unaveraged TestInput run by test1/3.
class TestOutput : public Step {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr, bool flag, bool clear,
             bool useComplement)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag),
        itsClear(clear),
        itsUseComplement(useComplement) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // The result of the selected and complement data depends on the
    // settings of itsFlag, itsClear, and itsUseComplement as follows:
    //    itsFlag itsUseComp itsClear     sel  comp
    //      T         T         T          T    F  *
    //      F         T         T          F    F
    //      T         F         T          F    T
    //      F         F         T          F    F
    //      T         T         F          T    T
    //      F         T         F          F    T  *
    //      T         F         F          T    T
    //      F         F         F          T    F
    // The lines marked with * are the cases in the if below.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    bool compFlag = itsFlag;
    bool selFlag = !itsClear;
    if (itsUseComplement && itsFlag == itsClear) {
      compFlag = !itsFlag;
      selFlag = itsClear;
    }
    result = compFlag;
    // All baselines except 2-2 should be selected.
    // Of them only channel 1,4,5 are selected.
    for (int i = 0; i < itsNBl; ++i) {
      if (i % 16 != 10) {
        for (int j = 0; j < itsNChan; ++j) {
          if (j == 1 || j == 4 || j == 5) {
            for (int k = 0; k < itsNCorr; ++k) {
              result(k, j, i) = selFlag;
            }
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
  bool itsFlag, itsClear, itsUseComplement;
};

// Test flagging a few antennae and freqs.
void test1(int ntime, int nbl, int nchan, int ncorr, bool flag, bool clear,
           bool useComplement) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  string mode;
  if (clear) {
    mode = (useComplement ? "clearComplement" : "clear");
  } else {
    mode = (useComplement ? "setComplement" : "set");
  }
  parset.add("mode", mode);
  parset.add("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add("baseline", "[rs01.*, *s*.*2, rs02.s01]");
  parset.add("countflag", "true");
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(new Counter(in, parset, "cnt"));
  Step::ShPtr step4(
      new TestOutput(ntime, nbl, nchan, ncorr, flag, clear, useComplement));
  dp3::steps::test::Execute({step1, step2, step3, step4});
}

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
    // A few baselines should be flagged (0, 7, 13, 15)
    // Furthermore channel 1,4,5,11,12,13 are flagged.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    result = false;
    for (int i = 0; i < itsNBl; ++i) {
      if (i % 16 == 0 || i % 16 == 7 || i % 16 == 13 || i % 16 == 15) {
        for (int j = 0; j < itsNChan; ++j) {
          if (j == 4) {
            for (int k = 0; k < itsNCorr; ++k) {
              result(k, j, i) = true;
            }
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

// Test flagging a few baselines, freqs, and channels.
void test2(int ntime, int nbl, int nchan, int ncorr) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, false);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add("chan", "[11..13, 4, 11, nchan/1000+1000..1000*nchan]");
  parset.add("baseline", "[[rs01.*,rs01.*],[*s*.*2,*s*.*2],[*s*.*2,rs02.*]]");
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test flagging a few antennae or freqs by using multiple steps.
void test3(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("expr", "s1");
  parset.add("s1.freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add("s1.expr", "s2");
  parset.add("s1.s2.baseline", "[rs01.*, *s*.*2, rs02.s01]");
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(
      new TestOutput(ntime, nbl, nchan, ncorr, flag, false, false));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Class to check result of flagged, unaveraged TestInput run by test4.
class TestOutput4 : public Step {
 public:
  TestOutput4(int ntime, int nbl, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // All baselines except autocorr should be flagged.
    // Furthermore channel 1,4,5 are flagged.
    casacore::Cube<bool> result(itsNCorr, itsNChan, itsNBl);
    result = true;
    for (int i = 0; i < itsNBl; ++i) {
      if (i % 16 == 0 || i % 16 == 5 || i % 16 == 10 || i % 16 == 15) {
        for (int j = 0; j < itsNChan; ++j) {
          if (j != 1 && j != 4 && j != 5) {
            for (int k = 0; k < itsNCorr; ++k) {
              result(k, j, i) = itsFlag;
            }
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
  bool itsFlag;
};

// Test flagging a few antennae and freqs by using multiple steps.
void test4(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("expr", "(s1&s1),(s2|s2)");
  parset.add("s1.freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add("s2.baseline", "[rs01.*, *s*.*2, rs02.s01]");
  parset.add("s2.corrtype", "cross");
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput4(ntime, nbl, nchan, ncorr, flag));
  dp3::steps::test::Execute({step1, step2, step3});
}

typedef bool CheckFunc(casacore::Complex value, double time, int ant1, int ant2,
                       const double* uvw);

// Class to check result of flagged, unaveraged TestInput run by test5.
class TestOutput5 : public Step {
 public:
  TestOutput5(CheckFunc* cfunc) : itsCount(0), itsCFunc(cfunc) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    const casacore::Cube<casacore::Complex>& data = buf.getData();
    const double* uvw = buf.getUVW().data();
    const casacore::IPosition& shp = data.shape();
    casacore::Cube<bool> result(shp);
    for (int i = 0; i < shp[2]; ++i) {
      int a1 = i / 4;
      int a2 = i % 4;
      for (int j = 0; j < shp[1]; ++j) {
        bool flag = false;
        for (int k = 0; k < shp[0]; ++k) {
          if (!flag)
            flag = itsCFunc(data(k, j, i), buf.getTime(), a1, a2, uvw + 3 * i);
        }
        for (int k = 0; k < shp[0]; ++k) {
          result(k, j, i) = flag;
        }
      }
    }
    BOOST_CHECK(allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {}

  int itsCount;
  CheckFunc* itsCFunc;
};

// Test flagging on a single parameter.
void test5(const string& key, const string& value, CheckFunc* cfunc) {
  // Create the steps.
  TestInput* in = new TestInput(2, 6, 5, 4, false);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add(key, value);
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput5(cfunc));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test flagging on multiple parameters.
void test6(const string& key1, const string& value1, const string& key2,
           const string& value2, CheckFunc* cfunc) {
  // Create the steps.
  TestInput* in = new TestInput(6, 10, 8, 4, false);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add(key1, value1);
  parset.add(key2, value2);
  Step::ShPtr step2(new PreFlagger(in, parset, ""));
  Step::ShPtr step3(new TestOutput5(cfunc));
  dp3::steps::test::Execute({step1, step2, step3});
}

bool checkBL(casacore::Complex, double, int a1, int a2, const double*) {
  return a1 == a2;
}
bool checkAmplMin(casacore::Complex data, double, int, int, const double*) {
  return abs(data) < 9.5;
}
bool checkAmplMax(casacore::Complex data, double, int, int, const double*) {
  return abs(data) > 31.5;
}
bool checkPhaseMin(casacore::Complex data, double, int, int, const double*) {
  return arg(data) < 1.4;
}
bool checkPhaseMax(casacore::Complex data, double, int, int, const double*) {
  return arg(data) > 2.1;
}
bool checkRealMin(casacore::Complex data, double, int, int, const double*) {
  return real(data) < 5.5;
}
bool checkRealMax(casacore::Complex data, double, int, int, const double*) {
  return real(data) > 29.4;
}
bool checkImagMin(casacore::Complex data, double, int, int, const double*) {
  return imag(data) < -1.4;
}
bool checkImagMax(casacore::Complex data, double, int, int, const double*) {
  return imag(data) > 20.5;
}
bool checkUVMin(casacore::Complex, double, int, int, const double* uvw) {
  return sqrt(uvw[0] * uvw[0] + uvw[1] * uvw[1]) <= 30;
}
bool checkUVBL(casacore::Complex, double, int a1, int a2, const double* uvw) {
  return sqrt(uvw[0] * uvw[0] + uvw[1] * uvw[1]) >= 30 && (a1 == 0 || a2 == 0);
}
bool checkBLMin(casacore::Complex, double, int a1, int a2, const double*) {
  return abs(a1 - a2) < 2;
}  // adjacent ant have bl<145
bool checkBLMinMax(casacore::Complex, double, int a1, int a2, const double*) {
  return abs(a1 - a2) != 1;
}  // adjacent ant have bl<145
bool checkTimeSlot(casacore::Complex, double time, int, int, const double*) {
  return time < 5;
}
bool checkNone(casacore::Complex, double, int, int, const double*) {
  return false;
}
bool checkAll(casacore::Complex, double, int, int, const double*) {
  return true;
}

// Test flagging on various fields.
void testMany() {
  test5("corrtype", "auto", &checkBL);
  test5("amplmin", "9.5", &checkAmplMin);
  test5("amplmax", "31.5", &checkAmplMax);
  test5("phasemin", "1.4", &checkPhaseMin);
  test5("phasemax", "2.1", &checkPhaseMax);
  test5("realmin", "5.5", &checkRealMin);
  test5("realmax", "29.4", &checkRealMax);
  test5("imagmin", "-1.4", &checkImagMin);
  test5("imagmax", "20.5", &checkImagMax);
  test5("uvmmin", "30", &checkUVMin);
  test6("uvmmax", "30", "baseline", "[rs01.s01]", &checkUVBL);
  test5("timeslot", "0", &checkTimeSlot);
  test5("abstime", "17-nov-1858/0:0:2..17nov1858/0:0:4", &checkTimeSlot);
  test5("abstime", "17-nov-1858/0:0:3+-1s", &checkTimeSlot);
  test5("reltime", "0:0:0..0:0:6", &checkTimeSlot);
  test5("reltime", "0:0:2+-1s", &checkTimeSlot);
  test5("reltime", "0:0:2+-20s", &checkAll);
  test6("abstime", "17-nov-1858/0:0:20+-19s", "reltime", "0:0:20+-20s",
        &checkAll);
  test6("abstime", "17-nov-1858/0:0:3+-2s", "reltime", "0:0:20+-20s",
        &checkTimeSlot);
  test6("abstime", "17-nov-1858/0:0:20+-19s", "reltime", "0:0:3+-2s",
        &checkTimeSlot);
  test6("abstime", "17-nov-1858/0:0:20+-9s", "reltime", "0:0:3+-2s",
        &checkNone);
  // Elevation is 12738s; azimuth=86121s
  test5("elevation", "180deg..190deg", &checkNone);
  test5("elevation", "12730s..12740s", &checkAll);
  test5("azimuth", "180deg..190deg", &checkNone);
  test5("azimuth", "86120s..86125s", &checkAll);
  test6("azimuth", "86120s..86125s", "elevation", "180deg..190deg", &checkNone);
  test6("azimuth", "86120s..86125s", "elevation", "12730s..12740s", &checkAll);
  test5("lst", "0.154d..0.155d", &checkAll);
  test5("blmin", "145", &checkBLMin);
  test6("blmin", "10", "blmax", "145", &checkBLMinMax);
}

BOOST_AUTO_TEST_CASE(test_1) { test1(10, 16, 32, 4, false, true, false); }

BOOST_AUTO_TEST_CASE(test_2) { test1(10, 16, 32, 4, true, true, false); }

BOOST_AUTO_TEST_CASE(test_3) { test1(10, 16, 32, 4, false, false, false); }

BOOST_AUTO_TEST_CASE(test_4) { test1(10, 16, 32, 4, true, false, false); }

BOOST_AUTO_TEST_CASE(test_5) { test1(10, 16, 32, 4, false, true, true); }

BOOST_AUTO_TEST_CASE(test_6) { test1(10, 16, 32, 4, true, true, true); }

BOOST_AUTO_TEST_CASE(test_7) { test1(10, 16, 32, 4, false, false, true); }

BOOST_AUTO_TEST_CASE(test_8) { test1(10, 16, 32, 4, true, false, true); }

BOOST_AUTO_TEST_CASE(test_9) { test2(2, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(test_10) { test2(2, 36, 16, 2); }

BOOST_AUTO_TEST_CASE(test_11) { test3(3, 16, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test_12) { test3(4, 16, 4, 2, true); }

BOOST_AUTO_TEST_CASE(test_13) { test4(3, 16, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test_many) { testMany(); }

BOOST_AUTO_TEST_SUITE_END()
