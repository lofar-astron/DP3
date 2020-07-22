// tPSet.cc: Test program for class PreFlagger::PSet
// Copyright (C) 2010
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$
//
// @author Ger van Diepen

#include "../../PreFlagger.h"
#include "../../../Common/ParameterSet.h"

#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Quanta/MVTime.h>

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

using DP3::DPPP::PreFlagger;

namespace {

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public DP3::DPPP::DPInput {
 public:
  TestInput(int nbl, int nchan, int ncorr) : itsNChan(nchan), itsNCorr(ncorr) {
    info().init(itsNCorr, 0, itsNChan, 0, 0, 50, string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(nbl);
    std::vector<int> ant2(nbl);
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
    std::vector<string> antNames{"rs01.s01", "rs02.s01", "cs01.s01",
                                 "cs01.s02"};
    std::vector<casacore::MPosition> antPos(4);
    std::vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    std::vector<double> chanWidth(nchan, 100000);
    casacore::Vector<double> chanFreqs(nchan);
    indgen(chanFreqs, 1050000., 100000.);
    info().set(chanFreqs, chanWidth);
  }

 private:
  virtual bool process(const DP3::DPPP::DPBuffer&) { return false; }
  virtual void finish() {}
  virtual void show(std::ostream&) const {}

  int itsNChan, itsNCorr;
};

}  // namespace

namespace DP3 {
namespace DPPP {
// This class name should match the friend in PreFlagger.
class TestPSet {
 public:
  TestPSet();

  void testNone();

  void testBL1();
  void testBL2();
  void testBL3();
  void testBLError1();
  void testBLError2();
  void testBLError3();

  void testChan1();
  void testChan2();
  void testChan3();

  void testTime1();
  void testTime2();
  void testTime3();
  void testTime4();
  void testTime5();
  void testTime6();
  void testTime7();

  void testMinMax1();
  void testMinMax2();
  void testMinMax3();
  void testMinMax4();
  void testMinMax5();

 private:
  std::unique_ptr<TestInput> in_;
  DP3::ParameterSet parset_;
};

TestPSet::TestPSet()
    : in_(boost::make_unique<TestInput>(16, 8, 4)), parset_() {}

void TestPSet::testNone() {
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(!(pset.itsFlagOnBL || pset.itsFlagOnAmpl || pset.itsFlagOnPhase ||
                pset.itsFlagOnReal || pset.itsFlagOnImag ||
                pset.itsFlagOnAzEl || pset.itsFlagOnUV));
}

void TestPSet::testBL1() {
  parset_.add("baseline", "[rs01.*, rs02.s01]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(!(pset.itsFlagOnAmpl || pset.itsFlagOnPhase ||
                pset.itsFlagOnReal || pset.itsFlagOnImag ||
                pset.itsFlagOnAzEl || pset.itsFlagOnUV) &&
              pset.itsFlagOnBL);
  // Make sure the matrix is correct.
  const casacore::Matrix<bool>& mat = pset.itsFlagBL;
  BOOST_CHECK_EQUAL(mat.shape(), casacore::IPosition(2, 4, 4));
  BOOST_CHECK(mat(0, 0) && mat(0, 1) && mat(0, 2) && mat(0, 3));
  BOOST_CHECK(mat(1, 0) && mat(1, 1) && mat(1, 2) && mat(1, 3));
  BOOST_CHECK(mat(2, 0) && mat(2, 1) && !mat(2, 2) && !mat(2, 3));
  BOOST_CHECK(mat(3, 0) && mat(3, 1) && !mat(3, 2) && !mat(3, 3));
}

void TestPSet::testBL2() {
  parset_.add("corrtype", "auto");
  parset_.add("baseline", "[rs01.*, [*s*.*2], rs02.s01]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  // Make sure the matrix is correct.
  const casacore::Matrix<bool>& mat = pset.itsFlagBL;
  BOOST_CHECK_EQUAL(mat.shape(), casacore::IPosition(2, 4, 4));
  BOOST_CHECK(mat(0, 0) && !mat(0, 1) && !mat(0, 2) && !mat(0, 3));
  BOOST_CHECK(!mat(1, 0) && mat(1, 1) && !mat(1, 2) && !mat(1, 3));
  BOOST_CHECK(!mat(2, 0) && !mat(2, 1) && !mat(2, 2) && !mat(2, 3));
  BOOST_CHECK(!mat(3, 0) && !mat(3, 1) && !mat(3, 2) && mat(3, 3));
}

void TestPSet::testBL3() {
  parset_.add("corrtype", "CROSS");
  parset_.add("baseline", "[[rs*, *s*.*1], [cs01.s01,cs01.s02]]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  // Make sure the matrix is correct.
  const casacore::Matrix<bool>& mat = pset.itsFlagBL;
  BOOST_CHECK_EQUAL(mat.shape(), casacore::IPosition(2, 4, 4));
  BOOST_CHECK(!mat(0, 0) && mat(0, 1) && mat(0, 2) && !mat(0, 3));
  BOOST_CHECK(mat(1, 0) && !mat(1, 1) && mat(1, 2) && !mat(1, 3));
  BOOST_CHECK(mat(2, 0) && mat(2, 1) && !mat(2, 2) && mat(2, 3));
  BOOST_CHECK(!mat(3, 0) && !mat(3, 1) && mat(3, 2) && !mat(3, 3));
}

void TestPSet::testBLError1() {
  parset_.add("corrtype", "crossx");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testBLError2() {
  parset_.add("baseline", "[[a,b,c]]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testBLError3() {
  parset_.add("baseline", "[[a,b], [ ] ]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testChan1() {
  in_ = boost::make_unique<TestInput>(16, 32, 4);
  parset_.add("chan", "[11..13, 4]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsChannels.size(), size_t{4});
  BOOST_CHECK_EQUAL(pset.itsChannels[0], size_t{4});
  BOOST_CHECK_EQUAL(pset.itsChannels[1], size_t{11});
  BOOST_CHECK_EQUAL(pset.itsChannels[2], size_t{12});
  BOOST_CHECK_EQUAL(pset.itsChannels[3], size_t{13});
  BOOST_CHECK_EQUAL(pset.itsChanFlags.shape(), casacore::IPosition(2, 4, 32));
  for (unsigned int i = 0; i < 32; ++i) {
    if (i == 4 || i == 11 || i == 12 || i == 13) {
      BOOST_CHECK(allEQ(pset.itsChanFlags.column(i), true));
    } else {
      BOOST_CHECK(allEQ(pset.itsChanFlags.column(i), false));
    }
  }
}

void TestPSet::testChan2() {
  in_ = boost::make_unique<TestInput>(16, 32, 4);
  parset_.add("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsChannels.size(), size_t{3});
  BOOST_CHECK_EQUAL(pset.itsChannels[0], size_t{1});
  BOOST_CHECK_EQUAL(pset.itsChannels[1], size_t{4});
  BOOST_CHECK_EQUAL(pset.itsChannels[2], size_t{5});
}

void TestPSet::testChan3() {
  in_ = boost::make_unique<TestInput>(16, 32, 4);
  parset_.add("chan", "[11..13, 4]");
  parset_.add("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsChannels.size(), size_t{1});
  BOOST_CHECK_EQUAL(pset.itsChannels[0], size_t{4});
}

void TestPSet::testTime1() {
  parset_.add("abstime", "[1mar2009/12:00:00..2mar2009/13:00:00]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsATimes.size(), size_t{2});
  casacore::Quantity q;
  casacore::MVTime::read(q, "1mar2009/12:00:00");
  BOOST_CHECK_EQUAL(q.getValue("s"), pset.itsATimes[0]);
  BOOST_CHECK_EQUAL(pset.itsATimes[1] - pset.itsATimes[0], 86400 + 3600);
}

void TestPSet::testTime2() {
  parset_.add("reltime", "[12:00:00..13:00:00, 16:00 +- 2min ]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsRTimes.size(), size_t{4});
  BOOST_CHECK_EQUAL(pset.itsRTimes[0], 12 * 3600);
  BOOST_CHECK_EQUAL(pset.itsRTimes[1], 13 * 3600);
  BOOST_CHECK_EQUAL(pset.itsRTimes[2], 16 * 3600 - 120);
  BOOST_CHECK_EQUAL(pset.itsRTimes[3], 16 * 3600 + 120);
}

void TestPSet::testTime3() {
  parset_.add("timeofday", "[22:00:00..2:00:00, 23:30 +- 1h ]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsTimes.size(), size_t{8});
  BOOST_CHECK_EQUAL(pset.itsTimes[0], -1);
  BOOST_CHECK_EQUAL(pset.itsTimes[1], 2 * 3600);
  BOOST_CHECK_EQUAL(pset.itsTimes[2], 22 * 3600);
  BOOST_CHECK_EQUAL(pset.itsTimes[3], 24 * 3600 + 1);
  BOOST_CHECK_EQUAL(pset.itsTimes[4], -1);
  BOOST_CHECK_EQUAL(pset.itsTimes[5], 1800);
  BOOST_CHECK_EQUAL(pset.itsTimes[6], 22 * 3600 + 1800);
  BOOST_CHECK_EQUAL(pset.itsTimes[7], 24 * 3600 + 1);
}

void TestPSet::testTime4() {
  parset_.add("timeslot", "[2..4, 10]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK_EQUAL(pset.itsTimeSlot.size(), size_t{4});
  BOOST_CHECK_EQUAL(pset.itsTimeSlot[0], 2u);
  BOOST_CHECK_EQUAL(pset.itsTimeSlot[1], 3u);
  BOOST_CHECK_EQUAL(pset.itsTimeSlot[2], 4u);
  BOOST_CHECK_EQUAL(pset.itsTimeSlot[3], 10u);
}

void TestPSet::testTime5() {
  parset_.add("reltime", "[12:00:00]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testTime6() {
  parset_.add("reltime", "[12:00:00..11:00:00]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testTime7() {
  parset_.add("abstime", "[12:00:00..13:00:00]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  BOOST_CHECK_THROW(pset.updateInfo(in_->getInfo()), std::exception);
}

void TestPSet::testMinMax1() {
  parset_.add("amplmin", "[23,,,45]");
  parset_.add("amplmax", "112.5");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(pset.itsFlagOnAmpl);
  BOOST_CHECK_EQUAL(pset.itsAmplMin.size(), size_t{4});
  BOOST_CHECK_EQUAL(pset.itsAmplMax.size(), size_t{4});
  BOOST_CHECK(casacore::near(pset.itsAmplMin[0], 23.));
  BOOST_CHECK(casacore::near(pset.itsAmplMin[1], -1e30));
  BOOST_CHECK(casacore::near(pset.itsAmplMin[2], -1e30));
  BOOST_CHECK(casacore::near(pset.itsAmplMin[3], 45.));
  BOOST_CHECK(casacore::near(pset.itsAmplMax[0], 112.5));
  BOOST_CHECK(casacore::near(pset.itsAmplMax[1], 112.5));
  BOOST_CHECK(casacore::near(pset.itsAmplMax[2], 112.5));
  BOOST_CHECK(casacore::near(pset.itsAmplMax[3], 112.5));
}

void TestPSet::testMinMax2() {
  parset_.add("phasemin", "[23]");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(pset.itsFlagOnPhase);
  BOOST_CHECK_EQUAL(pset.itsAmplMin.size(), size_t{4});
  BOOST_CHECK_EQUAL(pset.itsAmplMax.size(), size_t{4});
  BOOST_CHECK(casacore::near(pset.itsPhaseMin[0], 23.));
  BOOST_CHECK(casacore::near(pset.itsPhaseMin[1], -1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMin[2], -1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMin[3], -1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMax[0], 1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMax[1], 1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMax[2], 1e30));
  BOOST_CHECK(casacore::near(pset.itsPhaseMax[3], 1e30));
}

void TestPSet::testMinMax3() {
  parset_.add("uvmmin", "23");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(pset.itsFlagOnUV);
  BOOST_CHECK(casacore::near(pset.itsMinUV, 23. * 23.));
  BOOST_CHECK(casacore::near(pset.itsMaxUV, 1e30));
}

void TestPSet::testMinMax4() {
  parset_.add("uvmmax", "23");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(pset.itsFlagOnUV);
  BOOST_CHECK(pset.itsMinUV < 0.);
  BOOST_CHECK(casacore::near(pset.itsMaxUV, 23. * 23.));
}

void TestPSet::testMinMax5() {
  parset_.add("uvmmin", "23");
  parset_.add("uvmmax", "123");
  PreFlagger::PSet pset(in_.get(), parset_, "");
  pset.updateInfo(in_->getInfo());
  BOOST_CHECK(pset.itsFlagOnUV);
  BOOST_CHECK(casacore::near(pset.itsMinUV, 23. * 23.));
  BOOST_CHECK(casacore::near(pset.itsMaxUV, 123. * 123.));
}
}  // end namespace DPPP
}  // end namespace DP3

BOOST_AUTO_TEST_SUITE(pset)

using DP3::DPPP::TestPSet;

BOOST_AUTO_TEST_CASE(test_none) { TestPSet().testNone(); }

BOOST_AUTO_TEST_CASE(test_bl1) { TestPSet().testBL1(); }
BOOST_AUTO_TEST_CASE(test_bl2) { TestPSet().testBL2(); }
BOOST_AUTO_TEST_CASE(test_bl3) { TestPSet().testBL3(); }
BOOST_AUTO_TEST_CASE(test_bl_error1) { TestPSet().testBLError1(); }
BOOST_AUTO_TEST_CASE(test_bl_error2) { TestPSet().testBLError2(); }
BOOST_AUTO_TEST_CASE(test_bl_error3) { TestPSet().testBLError3(); }

BOOST_AUTO_TEST_CASE(test_chan1) { TestPSet().testChan1(); }
BOOST_AUTO_TEST_CASE(test_chan2) { TestPSet().testChan2(); }
BOOST_AUTO_TEST_CASE(test_chan3) { TestPSet().testChan3(); }

BOOST_AUTO_TEST_CASE(test_time1) { TestPSet().testTime1(); }
BOOST_AUTO_TEST_CASE(test_time2) { TestPSet().testTime2(); }
BOOST_AUTO_TEST_CASE(test_time3) { TestPSet().testTime3(); }
BOOST_AUTO_TEST_CASE(test_time4) { TestPSet().testTime4(); }
BOOST_AUTO_TEST_CASE(test_time5) { TestPSet().testTime5(); }
BOOST_AUTO_TEST_CASE(test_time6) { TestPSet().testTime6(); }
BOOST_AUTO_TEST_CASE(test_time7) { TestPSet().testTime7(); }

BOOST_AUTO_TEST_CASE(test_min_max1) { TestPSet().testMinMax1(); }
BOOST_AUTO_TEST_CASE(test_min_max2) { TestPSet().testMinMax2(); }
BOOST_AUTO_TEST_CASE(test_min_max3) { TestPSet().testMinMax3(); }
BOOST_AUTO_TEST_CASE(test_min_max4) { TestPSet().testMinMax4(); }
BOOST_AUTO_TEST_CASE(test_min_max5) { TestPSet().testMinMax5(); }

BOOST_AUTO_TEST_SUITE_END()
