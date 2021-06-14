// tTestDynStep.cc: Test program for class TestDynStep
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

// This file is mostly a copy of DPPP/test/tAverager.cc.
// Only the way the average steps are created is different.
// Average is used underneath TestDynStep.

#include "tStepCommon.h"
#include "../../base/DPBuffer.h"
#include "../../base/DPInfo.h"
#include "../../base/InputStep.h"
#include "../../base/DPRun.h"

#include "../../common/ParameterSet.h"
#include "../../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

using namespace LOFAR;
using namespace dp3::base;
using namespace casacore;

BOOST_AUTO_TEST_SUITE(dyn_step)

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
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
    }
    DPBuffer buf;
    buf.setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags(flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags(fullResFlags);
    Matrix<double> uvw(3, itsNBl);
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
    info().init(itsNCorr, 0, itsNChan, itsNTime, 100, 5, std::string(),
                std::string());
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
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    Cube<bool> fullResFlags(itsNChan, itsNAvgTime, itsNBl);
    fullResFlags = true;  // takes care of missing times at the end
    weights = 0;
    for (int j = itsCount * itsNAvgTime; j < itsCount * itsNAvgTime + navgtime;
         ++j) {
      for (int i = 0; i < int(data.size()); ++i) {
        data.data()[i] += Complex(i + j * 10, i - 1000 + j * 6);
        weights.data()[i] += float(1);
      }
      fullResFlags(Slicer(IPosition(3, 0, 0, 0),
                          IPosition(3, itsNChan, navgtime, itsNBl))) = itsFlag;
    }
    Cube<Complex> result(itsNCorr, nchan, itsNBl);
    Cube<float> resultw(itsNCorr, nchan, itsNBl);
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
    /// cout << imag(buf.getData()) << endl<<imag(result);
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), itsFlag));
    BOOST_CHECK(near(buf.getTime(), 2 + 5 * (itsCount * itsNAvgTime +
                                             (itsNAvgTime - 1) / 2.)));
    BOOST_CHECK(allNear(buf.getWeights(), resultw, 1e-5));
    if (navgtime == itsNAvgTime) {
      Matrix<double> uvw(3, itsNBl);
      indgen(uvw, 100 * (itsCount * itsNAvgTime + 0.5 * (itsNAvgTime - 1)));
      BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    }
    // cout <<buf.getFullResFlags()<< fullResFlags;
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(int(info.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(info.nchan()), 1 + (itsNChan - 1) / itsNAvgChan);
    BOOST_CHECK_EQUAL(int(info.ntime()), 1 + (itsNTime - 1) / itsNAvgTime);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5 * itsNAvgTime);
    BOOST_CHECK_EQUAL(int(info.nchanAvg()), itsNAvgChan);
    BOOST_CHECK_EQUAL(int(info.ntimeAvg()), itsNAvgTime);
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
    Cube<Complex> data(itsNrCorr, itsNrChan, itsNrBl);
    Cube<float> weights(itsNrCorr, itsNrChan, itsNrBl);
    Cube<bool> flags(itsNrCorr, itsNrChan, itsNrBl);
    int i = 0;
    for (int ib = 0; ib < itsNrBl; ++ib) {
      for (int ic = 0; ic < itsNrChan; ++ic) {
        for (int ip = 0; ip < itsNrCorr; ++ip) {
          data.data()[i] = Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
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
    Vector<rownr_t> rownrs(1, itsCount);
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
    info().init(itsNrCorr, 0, itsNrChan, itsNrTime, 100, 5, std::string(),
                std::string());
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNChan, 100000.);
    for (int i = 0; i < itsNChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }
  int itsCount, itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
  Cube<bool> itsFullResFlags;
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
    Cube<Complex> result(itsNrCorr, 1, itsNrBl);
    Cube<float> weights(itsNrCorr, 1, itsNrBl);
    Cube<bool> flags(itsNrCorr, 1, itsNrBl);
    Cube<bool> fullResFlags(itsNrChan, itsNrTime, itsNrBl);
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
                  weight * Complex(i + it * 10, i - 1000 + it * 6);
              weights(ip, 0, ib) += weight;
              ///  cout << result(ip,0,ib)  << weight << endl;
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
    /// cout << real(buf.getData()) << endl<<real(result);
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), flags));
    BOOST_CHECK(near(buf.getTime(), 2. + 5 * (itsNrTime - 1) / 2.));
    BOOST_CHECK(allNear(buf.getWeights(), weights, 1e-5));
    Matrix<double> uvw(3, itsNrBl);
    indgen(uvw);
    BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    /// cout <<buf.getFullResFlags()<< fullResFlags;
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(int(info.origNChan()), itsNrChan);
    BOOST_CHECK_EQUAL(info.nchan(), 1);
    BOOST_CHECK_EQUAL(info.ntime(), 1);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5 * itsNrTime);
    BOOST_CHECK_EQUAL(int(info.nchanAvg()), itsNrChan);
    BOOST_CHECK_EQUAL(int(info.ntimeAvg()), itsNrTime);
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
        /// cout << "flagged " <<itsCount <<' '<<  i << endl;
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
    Cube<Complex> result(itsNrCorr, 1, itsNrBl);
    Cube<float> weights(itsNrCorr, 1, itsNrBl);
    Cube<bool> flags(itsNrCorr, 1, itsNrBl);
    Cube<bool> fullResFlags(itsNrChan, itsNrTime, itsNrBl);
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
            /// cout << "out4 flagged "<< tf<<' '<< i/itsNrCorr<<' ' <<ib<<'
            /// '<<ic/2 << endl;
            i += itsNrCorr;
          } else {
            for (int ip = 0; ip < itsNrCorr; ++ip) {
              if ((it + 2 * ib + 3 * ic) % 7 != 0) {
                float weight = (1 + (it + ib + ic) % 5) / 5.;
                result(ip, 0, ib) +=
                    weight * Complex(i + it * 10, i - 1000 + it * 6);
                weights(ip, 0, ib) += weight;
                ///  cout << result(ip,0,ib)  << weight << endl;
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
    /// cout << real(buf.getData()) << endl<<real(result);
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(allEQ(buf.getFlags(), flags));
    BOOST_CHECK(near(buf.getTime(), 2. + 5 * (itsNrTime - 1) / 2.));
    BOOST_CHECK(allNear(buf.getWeights(), weights, 1e-5));
    Matrix<double> uvw(3, itsNrBl);
    indgen(uvw);
    BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-5));
    /// cout <<buf.getFullResFlags()<< fullResFlags;
    BOOST_CHECK(allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK_EQUAL(int(info.origNChan()), itsNrChan);
    BOOST_CHECK_EQUAL(info.nchan(), 1);
    BOOST_CHECK_EQUAL(info.ntime(), 1);
    BOOST_CHECK_EQUAL(info.timeInterval(), 5 * itsNrTime);
    BOOST_CHECK_EQUAL(int(info.nchanAvg()), itsNrChan);
    BOOST_CHECK_EQUAL(int(info.ntimeAvg()), itsNrTime);
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr, itsStep;
};

// Test simple averaging without flagged points.
void test1(int ntime, int nbl, int nchan, int ncorr, int navgtime, int navgchan,
           bool flag) {
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " navgtime=" << navgtime
       << " navgchan=" << navgchan << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  dp3::common::ParameterSet parset;
  parset.add("freqstep", std::to_string(navgchan));
  parset.add("timestep", std::to_string(navgtime));
  Step::ShPtr step2 = DPRun::findStepCtor("TestDynDPPP")(in, parset, "");
  Step::ShPtr step3(
      new TestOutput(ntime, nbl, nchan, ncorr, navgtime, navgchan, flag));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Like test1, but the averaging is done in two steps.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " navgtime=2"
       << " navgchan=4" << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  dp3::common::ParameterSet parset1, parset2;
  parset1.add("freqstep", "4");
  parset2.add("timestep", "2");
  Step::ShPtr step2a = DPRun::findStepCtor("TestDynDPPP")(in, parset1, "");
  Step::ShPtr step2b = DPRun::findStepCtor("TestDynDPPP")(in, parset2, "");
  Step::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr, 2, 4, flag));
  dp3::steps::test::Execute({step1, step2a, step2b, step3});
}

// Do tests with weighting and some flagged points.
void test3(int nrbl, int nrcorr) {
  {
    cout << "test3: ntime=2 nrbl=" << nrbl << " nchan=2 ncorr=" << nrcorr
         << endl;
    cout << "  navgtime=2 navgchan=2" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(2, nrbl, 2, nrcorr);
    Step::ShPtr step1(in);
    dp3::common::ParameterSet parset1;
    parset1.add("freqstep", "2");
    parset1.add("timestep", "2");
    Step::ShPtr step2a = DPRun::findStepCtor("TestDynDPPP")(in, parset1, "");
    Step::ShPtr step3(new TestOutput3(2, nrbl, 2, nrcorr));
    dp3::steps::test::Execute({step1, step2a, step3});
  }
  {
    cout << "test3: ntime=4 nrbl=" << nrbl << " nchan=8 ncorr=" << nrcorr
         << endl;
    cout << "  [navgtime=2 navgchan=4], [navgtime=2 navgchan=2]" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(4, nrbl, 8, nrcorr);
    Step::ShPtr step1(in);
    dp3::common::ParameterSet parset1, parset2;
    parset1.add("freqstep", "4");
    parset1.add("timestep", "2");
    parset2.add("freqstep", "2");
    parset2.add("timestep", "2");
    Step::ShPtr step2a = DPRun::findStepCtor("TestDynDPPP")(in, parset1, "");
    Step::ShPtr step2b = DPRun::findStepCtor("TestDynDPPP")(in, parset2, "");
    Step::ShPtr step3(new TestOutput3(4, nrbl, 8, nrcorr));
    dp3::steps::test::Execute({step1, step2a step2b, step3});
  }
}

// Do tests with averaging and flagging steps to see if the flags are
// promoted to the FULLRES flags.
void test4(int nrbl, int nrcorr, int flagstep) {
  {
    cout << "test4: ntime=4 nrbl=" << nrbl << " nchan=8 ncorr=" << nrcorr
         << endl;
    cout << "  [navgtime=2 navgchan=2], [flagstep=" << flagstep
         << "] [navgtime=2 navgchan=4]" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(4, nrbl, 8, nrcorr);
    Step::ShPtr step1(in);
    dp3::common::ParameterSet parset1, parset2;
    parset1.add("freqstep", "2");
    parset1.add("timestep", "2");
    parset2.add("freqstep", "4");
    parset2.add("timestep", "2");
    Step::ShPtr step2a = DPRun::findStepCtor("TestDynDPPP")(in, parset1, "");
    Step::ShPtr step2b(new TestFlagger(flagstep));
    Step::ShPtr step2c = DPRun::findStepCtor("TestDynDPPP")(in, parset2, "");
    Step::ShPtr step3(new TestOutput4(4, nrbl, 8, nrcorr, flagstep));
    dp3::steps::test::Execute({step1, step2a, step2b, step2c});
  }
}

BOOST_AUTO_TEST_CASE(test1a) { test1(10, 3, 32, 4, 2, 4, false); }

BOOST_AUTO_TEST_CASE(test1b) { test1(10, 3, 30, 1, 3, 3, true); }

BOOST_AUTO_TEST_CASE(test1c) { test1(10, 3, 30, 1, 3, 3, false); }

BOOST_AUTO_TEST_CASE(test1d) { test1(11, 3, 30, 2, 3, 3, false); }

BOOST_AUTO_TEST_CASE(test1e) { test1(10, 3, 32, 4, 1, 32, false); }

BOOST_AUTO_TEST_CASE(test1f) { test1(10, 3, 32, 1, 1, 1, false); }

BOOST_AUTO_TEST_CASE(test2a) { test2(10, 3, 32, 2, true); }

BOOST_AUTO_TEST_CASE(test2b) { test2(10, 3, 32, 2, false); }

BOOST_AUTO_TEST_CASE(test3a) { test3(1, 1); }

BOOST_AUTO_TEST_CASE(test3b) { test3(10, 4); }

BOOST_AUTO_TEST_CASE(test4a) { test4(1, 4, 3); }

BOOST_AUTO_TEST_CASE(test4b) { test4(20, 4, 5); }

BOOST_AUTO_TEST_SUITE_END()
