// tAverager.cc: Test program for class Averager
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../Averager.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Averager;
using dp3::steps::InputStep;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(averager)

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
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
          casacore::Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
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
    casacore::Matrix<double> uvw(3, itsNBl);
    indgen(uvw, double(itsCount * 100));
    buf.setUVW(uvw);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {
    // Use timeInterval=5
    info().init(itsNCorr, 0, itsNChan, itsNTime, 100, 5, string(), string());
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNChan, 100000.);
    for (int i = 0; i < itsNChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }
  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of averaging TestInput.
class TestOutput : public Step {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr, int navgtime,
             int navgchan, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsNAvgTime(navgtime),
        itsNAvgChan(navgchan),
        itsFlag(flag) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    int nchan = 1 + (itsNChan - 1) / itsNAvgChan;
    int navgtime = std::min(itsNAvgTime, itsNTime - itsCount * itsNAvgTime);
    // Fill expected result in similar way as TestInput.
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    casacore::Cube<bool> fullResFlags(itsNChan, itsNAvgTime, itsNBl);
    fullResFlags = true;  // takes care of missing times at the end
    weights = 0;
    for (int j = itsCount * itsNAvgTime; j < itsCount * itsNAvgTime + navgtime;
         ++j) {
      for (int i = 0; i < int(data.size()); ++i) {
        data.data()[i] += casacore::Complex(i + j * 10, i - 1000 + j * 6);
        weights.data()[i] += float(1);
      }
      fullResFlags(casacore::Slicer(
          casacore::IPosition(3, 0, 0, 0),
          casacore::IPosition(3, itsNChan, navgtime, itsNBl))) = itsFlag;
    }
    casacore::Cube<casacore::Complex> result(itsNCorr, nchan, itsNBl);
    casacore::Cube<float> resultw(itsNCorr, nchan, itsNBl);
    resultw = 0;
    // Average to get the true expected result.
    for (int k = 0; k < itsNBl; ++k) {
      for (int i = 0; i < itsNCorr; ++i) {
        for (int j = 0; j < nchan; ++j) {
          int jc;
          for (jc = j * itsNAvgChan;
               jc < std::min((j + 1) * itsNAvgChan, itsNChan); ++jc) {
            result(i, j, k) += data(i, jc, k);
            resultw(i, j, k) += weights(i, jc, k);
          }
          result(i, j, k) /= float(navgtime * (jc - j * itsNAvgChan));
        }
      }
    }
    // Check the averaged result.
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), itsFlag));
    BOOST_CHECK(casacore::near(
        buf.getTime(),
        2 + 5 * (itsCount * itsNAvgTime + (itsNAvgTime - 1) / 2.)));
    BOOST_CHECK(allNear(buf.getWeights(), resultw, 1e-5));
    if (navgtime == itsNAvgTime) {
      casacore::Matrix<double> uvw(3, itsNBl);
      indgen(uvw, 100 * (itsCount * itsNAvgTime + 0.5 * (itsNAvgTime - 1)));
      BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    }
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(itsNChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL(1 + (itsNChan - 1) / itsNAvgChan, int(info.nchan()));
    BOOST_CHECK_EQUAL(1 + (itsNTime - 1) / itsNAvgTime, int(info.ntime()));
    BOOST_CHECK_EQUAL(5 * itsNAvgTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNAvgChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL(itsNAvgTime, int(info.ntimeAvg()));
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
  bool itsFlag;
};

// More elaborate class which can set different flags and weights.
class TestInput3 : public InputStep {
 public:
  TestInput3(int nrtime, int nrbl, int nrchan, int nrcorr)
      : itsCount(0),
        itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr) {
    itsFullResFlags.resize(itsNrChan, 1, nrbl);
  }

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (itsCount == itsNrTime) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(itsNrCorr, itsNrChan, itsNrBl);
    casacore::Cube<float> weights(itsNrCorr, itsNrChan, itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr, itsNrChan, itsNrBl);
    int i = 0;
    for (int ib = 0; ib < itsNrBl; ++ib) {
      for (int ic = 0; ic < itsNrChan; ++ic) {
        for (int ip = 0; ip < itsNrCorr; ++ip) {
          data.data()[i] =
              casacore::Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
          weights.data()[i] = (1 + (itsCount + ib + ic) % 5) / 5.;
          flags.data()[i] = ((itsCount + 2 * ib + 3 * ic) % 7 == 0);
          i++;
        }
        itsFullResFlags(ic, 0, ib) = ((itsCount + 2 * ib + 3 * ic) % 7 == 0);
      }
    }
    DPBuffer buf;
    buf.setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    buf.setWeights(weights);
    buf.setFlags(flags);
    vector<dp3::common::rownr_t> rownrs(1, itsCount);
    buf.setRowNrs(rownrs);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void getUVW(const casacore::RefRows&, double, DPBuffer& buf) {
    buf.getUVW().resize(3, itsNrBl);
    indgen(buf.getUVW());
  }
  virtual bool getFullResFlags(const casacore::RefRows&, DPBuffer& buf) {
    buf.getFullResFlags().assign(itsFullResFlags);
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {
    // Use timeInterval=5
    info().init(itsNrCorr, 0, itsNrChan, itsNrTime, 100, 5, string(), string());
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNrChan, 100000.);
    for (int i = 0; i < itsNrChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }
  int itsCount, itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
  casacore::Cube<bool> itsFullResFlags;
};

// Class to check result of averaging TestInput3.
// All input must be averaged (in one or more steps) to a single value
// per corr/baseline.
class TestOutput3 : public Step {
 public:
  TestOutput3(int nrtime, int nrbl, int nrchan, int nrcorr)
      : itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    casacore::Cube<casacore::Complex> result(itsNrCorr, 1, itsNrBl);
    casacore::Cube<float> weights(itsNrCorr, 1, itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr, 1, itsNrBl);
    casacore::Cube<bool> fullResFlags(itsNrChan, itsNrTime, itsNrBl);
    weights = float(0);
    flags = true;
    fullResFlags = true;
    // Create data in the same way as in TestInput3.
    for (int it = 0; it < itsNrTime; ++it) {
      int i = 0;
      for (int ib = 0; ib < itsNrBl; ++ib) {
        for (int ic = 0; ic < itsNrChan; ++ic) {
          for (int ip = 0; ip < itsNrCorr; ++ip) {
            if ((it + 2 * ib + 3 * ic) % 7 != 0) {
              float weight = (1 + (it + ib + ic) % 5) / 5.;
              result(ip, 0, ib) +=
                  weight * casacore::Complex(i + it * 10, i - 1000 + it * 6);
              weights(ip, 0, ib) += weight;
              flags(ip, 0, ib) = false;
              fullResFlags(ic, it, ib) = false;
            }
            i++;
          }
        }
      }
    }
    BOOST_CHECK(allNE(weights, float(0.)));
    for (unsigned int i = 0; i < result.size(); ++i) {
      result.data()[i] /= weights.data()[i];
    }
    // Check the averaged result.
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), flags));
    BOOST_CHECK(casacore::near(buf.getTime(), 2. + 5 * (itsNrTime - 1) / 2.));
    BOOST_CHECK(allNear(buf.getWeights(), weights, 1e-5));
    casacore::Matrix<double> uvw(3, itsNrBl);
    indgen(uvw);
    BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(itsNrChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL(size_t{1}, info.nchan());
    BOOST_CHECK_EQUAL(size_t{1}, info.ntime());
    BOOST_CHECK_EQUAL(5 * itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNrChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL(itsNrTime, int(info.ntimeAvg()));
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
};

// Simple class to flag every step-th XX point.
class TestFlagger : public Step {
 public:
  TestFlagger(int step) : itsCount(0), itsStep(step) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    DPBuffer buf2(buf);
    int ncorr = buf2.getFlags().shape()[0];
    int np = buf2.getFlags().size() / ncorr;
    bool* flagPtr = buf2.getFlags().data();
    for (int i = 0; i < np; ++i) {
      if ((i + itsCount) % itsStep == 0) {
        for (int j = 0; j < ncorr; ++j) {
          flagPtr[i * ncorr + j] = true;
        }
      }
    }
    getNextStep()->process(buf2);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}

  int itsCount, itsStep;
};

// Class to check result of averaging and flagging TestInput3.
// First the data are averaged from 8,4 to 4,2, then every step-th point
// is flagged, and finally it is averaged to 1,1.
class TestOutput4 : public Step {
 public:
  TestOutput4(int nrtime, int nrbl, int nrchan, int nrcorr, int step)
      : itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr),
        itsStep(step) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    casacore::Cube<casacore::Complex> result(itsNrCorr, 1, itsNrBl);
    casacore::Cube<float> weights(itsNrCorr, 1, itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr, 1, itsNrBl);
    casacore::Cube<bool> fullResFlags(itsNrChan, itsNrTime, itsNrBl);
    weights = float(0);
    flags = true;
    fullResFlags = true;
    // Create data in the same way as in TestInput3.
    for (int it = 0; it < itsNrTime; ++it) {
      int i = 0;
      for (int ib = 0; ib < itsNrBl; ++ib) {
        for (int ic = 0; ic < itsNrChan; ++ic) {
          // TestFlagger flags every step-th point of 2x2 averaged data.
          int tf = it / 2;  // same as itsCount in testFlagger
          if (((ib * itsNrChan + ic) / 2 + tf) % itsStep == 0) {
            i += itsNrCorr;
          } else {
            for (int ip = 0; ip < itsNrCorr; ++ip) {
              if ((it + 2 * ib + 3 * ic) % 7 != 0) {
                float weight = (1 + (it + ib + ic) % 5) / 5.;
                result(ip, 0, ib) +=
                    weight * casacore::Complex(i + it * 10, i - 1000 + it * 6);
                weights(ip, 0, ib) += weight;
                flags(ip, 0, ib) = false;
                fullResFlags(ic, it, ib) = false;
              }
              i++;
            }
          }
        }
      }
    }
    for (unsigned int i = 0; i < result.size(); ++i) {
      if (!flags.data()[i]) {
        result.data()[i] /= weights.data()[i];
      }
    }
    // Check the averaged result.
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), flags));
    BOOST_CHECK(casacore::near(buf.getTime(), 2. + 5 * (itsNrTime - 1) / 2.));
    BOOST_CHECK(allNear(buf.getWeights(), weights, 1e-5));
    casacore::Matrix<double> uvw(3, itsNrBl);
    indgen(uvw);
    BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(itsNrChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL(size_t{1}, info.nchan());
    BOOST_CHECK_EQUAL(size_t{1}, info.ntime());
    BOOST_CHECK_EQUAL(5 * itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNrChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL(itsNrTime, int(info.ntimeAvg()));
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr, itsStep;
};

// Test simple averaging without flagged points.
void test1(int ntime, int nbl, int nchan, int ncorr, int navgtime, int navgchan,
           bool flag) {
  // Create the steps.
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqstep", std::to_string(navgchan));
  parset.add("timestep", std::to_string(navgtime));
  auto step2 = std::make_shared<Averager>(*step1, parset, "");
  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, navgtime,
                                            navgchan, flag);
  dp3::steps::test::Execute({step1, step2, step3});
}

// Like test 1, but specify target resolution
void test1resolution(int ntime, int nbl, int nchan, int ncorr,
                     double timeresolution, double freqresolution,
                     string frequnit, bool flag) {
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqresolution", std::to_string(freqresolution) + frequnit);
  parset.add("timeresolution", std::to_string(timeresolution));
  auto step2 = std::make_shared<Averager>(*step1, parset, "");

  if (!frequnit.empty()) {
    casacore::Quantity q(freqresolution, frequnit);
    freqresolution = q.getValue("Hz", true);
  }

  int navgchan = std::max(1, int(freqresolution / 100000 + 0.5));
  int navgtime = std::max(1, int(timeresolution / 5. + 0.5));

  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, navgtime,
                                            navgchan, flag);
  dp3::steps::test::Execute({step1, step2, step3});
}

// Like test1, but the averaging is done in two steps.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset1, parset2;
  parset1.add("freqstep", "4");
  parset2.add("timestep", "2");
  auto step2a = std::make_shared<Averager>(*step1, parset1, "");
  auto step2b = std::make_shared<Averager>(*step1, parset2, "");
  auto step3 =
      std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, 2, 4, flag);
  dp3::steps::test::Execute({step1, step2a, step2b, step3});
}

// Do tests with weighting and some flagged points.
void test3(int nrbl, int nrcorr) {
  {
    // Create the steps.
    auto step1 = std::make_shared<TestInput3>(2, nrbl, 2, nrcorr);
    ParameterSet parset1;
    parset1.add("freqstep", "2");
    parset1.add("timestep", "2");
    auto step2 = std::make_shared<Averager>(*step1, parset1, "");
    auto step3 = std::make_shared<TestOutput3>(2, nrbl, 2, nrcorr);
    dp3::steps::test::Execute({step1, step2, step3});
  }
  {
    // Create the steps.
    auto step1 = std::make_shared<TestInput3>(4, nrbl, 8, nrcorr);
    ParameterSet parset1, parset2;
    parset1.add("freqstep", "4");
    parset1.add("timestep", "2");
    parset2.add("freqstep", "2");
    parset2.add("timestep", "2");
    auto step2a = std::make_shared<Averager>(*step1, parset1, "");
    auto step2b = std::make_shared<Averager>(*step1, parset2, "");
    auto step3 = std::make_shared<TestOutput3>(4, nrbl, 8, nrcorr);
    dp3::steps::test::Execute({step1, step2a, step2b, step3});
  }
}

// Do tests with averaging and flagging steps to see if the flags are
// promoted to the FULLRES flags.
void test4(int nrbl, int nrcorr, int flagstep) {
  {
    // Create the steps.
    auto step1 = std::make_shared<TestInput3>(4, nrbl, 8, nrcorr);
    ParameterSet parset1, parset2;
    parset1.add("freqstep", "2");
    parset1.add("timestep", "2");
    parset2.add("freqstep", "4");
    parset2.add("timestep", "2");
    auto step2a = std::make_shared<Averager>(*step1, parset1, "");
    auto step2b = std::make_shared<TestFlagger>(flagstep);
    auto step2c = std::make_shared<Averager>(*step1, parset2, "");
    auto step3 = std::make_shared<TestOutput4>(4, nrbl, 8, nrcorr, flagstep);
    dp3::steps::test::Execute({step1, step2a, step2b, step2c, step3});
  }
}

BOOST_AUTO_TEST_CASE(testaverager1) { test1(10, 3, 32, 4, 2, 4, false); }

BOOST_AUTO_TEST_CASE(testaverager2) { test1(10, 3, 30, 1, 3, 3, true); }

BOOST_AUTO_TEST_CASE(testaverager3) { test1(10, 3, 30, 1, 3, 3, false); }
BOOST_AUTO_TEST_CASE(testaverager4) { test1(11, 3, 30, 2, 3, 3, false); }

BOOST_AUTO_TEST_CASE(testaverager5) { test1(10, 3, 32, 4, 1, 32, false); }

BOOST_AUTO_TEST_CASE(testaverager6) { test1(10, 3, 32, 1, 1, 1, false); }

BOOST_AUTO_TEST_CASE(testaverager7) { test2(10, 3, 32, 2, true); }

BOOST_AUTO_TEST_CASE(testaverager8) { test2(10, 3, 32, 2, false); }

BOOST_AUTO_TEST_CASE(testaverager9) { test3(1, 1); }

BOOST_AUTO_TEST_CASE(testaverager10) { test3(10, 4); }

BOOST_AUTO_TEST_CASE(testaverager11) { test4(1, 4, 3); }

BOOST_AUTO_TEST_CASE(testaverager12) { test4(20, 4, 5); }

BOOST_AUTO_TEST_CASE(testresolution1) {
  test1resolution(10, 3, 32, 4, 10., 100000, "Hz", false);
}

BOOST_AUTO_TEST_CASE(testresolution2) {
  test1resolution(11, 3, 32, 4, 1., 800, "kHz", false);
}

BOOST_AUTO_TEST_CASE(testresolution3) {
  test1resolution(11, 3, 32, 4, 15., 0.4, "MHz", false);
}

BOOST_AUTO_TEST_SUITE_END()
