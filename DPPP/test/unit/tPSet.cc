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
#include "../../DPInfo.h"
#include "../../../Common/ParameterSet.h"

#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Quanta/MVTime.h>

#include <boost/test/unit_test.hpp>

#include <iostream>

namespace {

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DP3::DPPP::DPInput
{
public:
  TestInput(int nbl, int nchan, int ncorr)
    : itsNChan(nchan), itsNCorr(ncorr)
  {
    info().init (itsNCorr, 0, itsNChan, 0, 0, 50, string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(nbl);
    std::vector<int> ant2(nbl);
    int st1 = 0;
    int st2 = 0;
    for (int i=0; i<nbl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<string> antNames {"rs01.s01", "rs02.s01", "cs01.s01", "cs01.s02"};
    std::vector<casacore::MPosition> antPos(4);
    std::vector<double> antDiam(4, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    std::vector<double> chanWidth(nchan, 100000);
    casacore::Vector<double> chanFreqs(nchan);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
private:
  virtual bool process (const DP3::DPPP::DPBuffer&) { return false; }
  virtual void finish() {}
  virtual void show (std::ostream&) const {}

  int itsNChan, itsNCorr;
};

} // end namespace anonymous

namespace DP3 {
  namespace DPPP {
    // This class name should match the friend in PreFlagger.
    class TestPSet
    {
    public:
      static void testNone();
      static void testBL();
      static void testChan();
      static void testTime();
      static void testMinMax();
    };

    void TestPSet::testNone()
     {
       TestInput* in = new TestInput(16, 8, 4);
       DPStep::ShPtr step1(in);
       ParameterSet parset;
       PreFlagger::PSet pset (in, parset, "");
       pset.updateInfo (in->getInfo());
       BOOST_CHECK (!(pset.itsFlagOnBL   || pset.itsFlagOnAmpl || pset.itsFlagOnPhase ||
                 pset.itsFlagOnReal || pset.itsFlagOnImag ||
                 pset.itsFlagOnAzEl || pset.itsFlagOnUV));
     }

    void TestPSet::testBL()
    {
      TestInput* in = new TestInput(16, 8, 4);
      DPStep::ShPtr step1(in);
      {
        ParameterSet parset;
        parset.add ("baseline", "[rs01.*, rs02.s01]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (!(pset.itsFlagOnAmpl || pset.itsFlagOnPhase || pset.itsFlagOnReal ||
                  pset.itsFlagOnImag || pset.itsFlagOnAzEl  || pset.itsFlagOnUV) &&
                pset.itsFlagOnBL);
        // Make sure the matrix is correct.
        const casacore::Matrix<bool>& mat = pset.itsFlagBL;
        BOOST_CHECK_EQUAL (mat.shape(), casacore::IPosition(2,4,4));
        BOOST_CHECK ( mat(0,0) &&  mat(0,1) &&  mat(0,2) &&  mat(0,3));
        BOOST_CHECK ( mat(1,0) &&  mat(1,1) &&  mat(1,2) &&  mat(1,3));
        BOOST_CHECK ( mat(2,0) &&  mat(2,1) && !mat(2,2) && !mat(2,3));
        BOOST_CHECK ( mat(3,0) &&  mat(3,1) && !mat(3,2) && !mat(3,3));
      }
      {
        ParameterSet parset;
        parset.add ("corrtype", "auto");
        parset.add ("baseline", "[rs01.*, [*s*.*2], rs02.s01]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        // Make sure the matrix is correct.
        const casacore::Matrix<bool>& mat = pset.itsFlagBL;
        BOOST_CHECK_EQUAL (mat.shape(), casacore::IPosition(2,4,4));
        BOOST_CHECK ( mat(0,0) && !mat(0,1) && !mat(0,2) && !mat(0,3));
        BOOST_CHECK (!mat(1,0) &&  mat(1,1) && !mat(1,2) && !mat(1,3));
        BOOST_CHECK (!mat(2,0) && !mat(2,1) && !mat(2,2) && !mat(2,3));
        BOOST_CHECK (!mat(3,0) && !mat(3,1) && !mat(3,2) &&  mat(3,3));
      }
      {
        ParameterSet parset;
        parset.add ("corrtype", "CROSS");
        parset.add ("baseline", "[[rs*, *s*.*1], [cs01.s01,cs01.s02]]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        // Make sure the matrix is correct.
        const casacore::Matrix<bool>& mat = pset.itsFlagBL;
        BOOST_CHECK_EQUAL (mat.shape(), casacore::IPosition(2,4,4));
        BOOST_CHECK (!mat(0,0) &&  mat(0,1) &&  mat(0,2) && !mat(0,3));
        BOOST_CHECK ( mat(1,0) && !mat(1,1) &&  mat(1,2) && !mat(1,3));
        BOOST_CHECK ( mat(2,0) &&  mat(2,1) && !mat(2,2) &&  mat(2,3));
        BOOST_CHECK (!mat(3,0) && !mat(3,1) &&  mat(3,2) && !mat(3,3));
      }
      // Some erronous ones.
      bool err = false;
      try {
        ParameterSet parset;
        parset.add ("corrtype", "crossx");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
      err = false;
      try {
        ParameterSet parset;
        parset.add ("baseline", "[[a,b,c]]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
      err = false;
      try {
        ParameterSet parset;
        parset.add ("baseline", "[[a,b], [ ] ]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
    }

    void TestPSet::testChan()
    {
      TestInput* in = new TestInput(16, 32, 4);
      DPStep::ShPtr step1(in);
      {
        ParameterSet parset;
        parset.add ("chan", "[11..13, 4]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsChannels.size(), size_t{4});
        BOOST_CHECK_EQUAL (pset.itsChannels[0], size_t{4});
        BOOST_CHECK_EQUAL (pset.itsChannels[1], size_t{11});
        BOOST_CHECK_EQUAL (pset.itsChannels[2], size_t{12});
        BOOST_CHECK_EQUAL (pset.itsChannels[3], size_t{13});
        BOOST_CHECK_EQUAL (pset.itsChanFlags.shape(), casacore::IPosition(2,4,32));
        for (unsigned int i=0; i<32; ++i) {
          if (i==4 || i==11 || i==12 || i==13) {
            BOOST_CHECK (allEQ(pset.itsChanFlags.column(i), true));
          } else {
            BOOST_CHECK (allEQ(pset.itsChanFlags.column(i), false));
          }
        }
      }
      {
        ParameterSet parset;
        parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsChannels.size(), size_t{3});
        BOOST_CHECK_EQUAL (pset.itsChannels[0], size_t{1});
        BOOST_CHECK_EQUAL (pset.itsChannels[1], size_t{4});
        BOOST_CHECK_EQUAL (pset.itsChannels[2], size_t{5});
      }
      {
        ParameterSet parset;
        parset.add ("chan", "[11..13, 4]");
        parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsChannels.size(), size_t{1});
        BOOST_CHECK_EQUAL (pset.itsChannels[0], size_t{4});
      }
    }

    void TestPSet::testTime()
    {
      TestInput* in = new TestInput(16, 8, 4);
      DPStep::ShPtr step1(in);
      {
        ParameterSet parset;
        parset.add ("abstime", "[1mar2009/12:00:00..2mar2009/13:00:00]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsATimes.size(), size_t{2});
        casacore::Quantity q;
        casacore::MVTime::read (q, "1mar2009/12:00:00");
        BOOST_CHECK_EQUAL (q.getValue("s"), pset.itsATimes[0]);
        BOOST_CHECK_EQUAL (pset.itsATimes[1] - pset.itsATimes[0], 86400+3600);
      }
      {
        ParameterSet parset;
        parset.add ("reltime", "[12:00:00..13:00:00, 16:00 +- 2min ]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsRTimes.size(), size_t{4});
        BOOST_CHECK_EQUAL (pset.itsRTimes[0], 12*3600);
        BOOST_CHECK_EQUAL (pset.itsRTimes[1], 13*3600);
        BOOST_CHECK_EQUAL (pset.itsRTimes[2], 16*3600-120);
        BOOST_CHECK_EQUAL (pset.itsRTimes[3], 16*3600+120);
      }
      {
        ParameterSet parset;
        parset.add ("timeofday", "[22:00:00..2:00:00, 23:30 +- 1h ]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsTimes.size(), size_t{8});
        BOOST_CHECK_EQUAL (pset.itsTimes[0], -1);
        BOOST_CHECK_EQUAL (pset.itsTimes[1], 2*3600);
        BOOST_CHECK_EQUAL (pset.itsTimes[2], 22*3600);
        BOOST_CHECK_EQUAL (pset.itsTimes[3], 24*3600+1);
        BOOST_CHECK_EQUAL (pset.itsTimes[4], -1);
        BOOST_CHECK_EQUAL (pset.itsTimes[5], 1800);
        BOOST_CHECK_EQUAL (pset.itsTimes[6], 22*3600+1800);
        BOOST_CHECK_EQUAL (pset.itsTimes[7], 24*3600+1);
      }
      {
        ParameterSet parset;
        parset.add ("timeslot", "[2..4, 10]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK_EQUAL (pset.itsTimeSlot.size(), size_t{4});
        BOOST_CHECK_EQUAL (pset.itsTimeSlot[0], 2u);
        BOOST_CHECK_EQUAL (pset.itsTimeSlot[1], 3u);
        BOOST_CHECK_EQUAL (pset.itsTimeSlot[2], 4u);
        BOOST_CHECK_EQUAL (pset.itsTimeSlot[3], 10u);
      }
      // Some erronous ones.
      bool err = false;
      try {
        ParameterSet parset;
        parset.add ("reltime", "[12:00:00]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
      err = false;
      try {
        ParameterSet parset;
        parset.add ("reltime", "[12:00:00..11:00:00]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
      err = false;
      try {
        ParameterSet parset;
        parset.add ("abstime", "[12:00:00..13:00:00]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
      } catch (std::exception& x) {
        err = true;
      }
      BOOST_CHECK (err);
    }

    void TestPSet::testMinMax()
    {
      using casacore::near;

      TestInput* in = new TestInput(16, 8, 4);
      DPStep::ShPtr step1(in);
      {
        ParameterSet parset;
        parset.add ("amplmin", "[23,,,45]");
        parset.add ("amplmax", "112.5");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (pset.itsFlagOnAmpl);
        BOOST_CHECK_EQUAL (pset.itsAmplMin.size(), size_t{4});
        BOOST_CHECK_EQUAL (pset.itsAmplMax.size(), size_t{4});
        BOOST_CHECK (near(pset.itsAmplMin[0], 23.));
        BOOST_CHECK (near(pset.itsAmplMin[1], -1e30));
        BOOST_CHECK (near(pset.itsAmplMin[2], -1e30));
        BOOST_CHECK (near(pset.itsAmplMin[3], 45.));
        BOOST_CHECK (near(pset.itsAmplMax[0], 112.5));
        BOOST_CHECK (near(pset.itsAmplMax[1], 112.5));
        BOOST_CHECK (near(pset.itsAmplMax[2], 112.5));
        BOOST_CHECK (near(pset.itsAmplMax[3], 112.5));
      }
      {
        ParameterSet parset;
        parset.add ("phasemin", "[23]");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (pset.itsFlagOnPhase);
        BOOST_CHECK_EQUAL (pset.itsAmplMin.size(), size_t{4});
        BOOST_CHECK_EQUAL (pset.itsAmplMax.size(), size_t{4});
        BOOST_CHECK (near(pset.itsPhaseMin[0], 23.));
        BOOST_CHECK (near(pset.itsPhaseMin[1], -1e30));
        BOOST_CHECK (near(pset.itsPhaseMin[2], -1e30));
        BOOST_CHECK (near(pset.itsPhaseMin[3], -1e30));
        BOOST_CHECK (near(pset.itsPhaseMax[0], 1e30));
        BOOST_CHECK (near(pset.itsPhaseMax[1], 1e30));
        BOOST_CHECK (near(pset.itsPhaseMax[2], 1e30));
        BOOST_CHECK (near(pset.itsPhaseMax[3], 1e30));
      }
      {
        ParameterSet parset;
        parset.add ("uvmmin", "23");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (pset.itsFlagOnUV);
        BOOST_CHECK (near(pset.itsMinUV, 23.*23.));
        BOOST_CHECK (near(pset.itsMaxUV, 1e30));
      }
      {
        ParameterSet parset;
        parset.add ("uvmmax", "23");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (pset.itsFlagOnUV);
        BOOST_CHECK (pset.itsMinUV < 0.);
        BOOST_CHECK (near(pset.itsMaxUV, 23.*23.));
      }
      {
        ParameterSet parset;
        parset.add ("uvmmin", "23");
        parset.add ("uvmmax", "123");
        PreFlagger::PSet pset (in, parset, "");
        pset.updateInfo (in->getInfo());
        BOOST_CHECK (pset.itsFlagOnUV);
        BOOST_CHECK (near(pset.itsMinUV, 23.*23.));
        BOOST_CHECK (near(pset.itsMaxUV, 123.*123.));
      }
    }
  } // end namespace DPPP
} // end namespace DP3

BOOST_AUTO_TEST_SUITE(pset)

using DP3::DPPP::TestPSet;

BOOST_AUTO_TEST_CASE( test_none ) {
  TestPSet::testNone();
}

BOOST_AUTO_TEST_CASE( test_bl ) {
  TestPSet::testBL();
}

BOOST_AUTO_TEST_CASE( test_chan ) {
  TestPSet::testChan();
}

BOOST_AUTO_TEST_CASE( test_time ) {
  TestPSet::testTime();
}

BOOST_AUTO_TEST_CASE( test_min_max ) {
  TestPSet::testMinMax();
}

BOOST_AUTO_TEST_SUITE_END()
