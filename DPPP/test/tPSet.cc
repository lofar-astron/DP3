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

#include <lofar_config.h>
#include <DPPP/PreFlagger.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <iostream>

using namespace DP3;
using namespace DP3::DPPP;
using namespace casacore;
using namespace std;

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int nbl, int nchan, int ncorr)
    : itsNChan(nchan), itsNCorr(ncorr)
  {
    info().init (itsNCorr, itsNChan, 0, 0, 50, string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    Vector<Int> ant1(nbl);
    Vector<Int> ant2(nbl);
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
    Vector<String> antNames(4);
    antNames[0] = "rs01.s01";
    antNames[1] = "rs02.s01";
    antNames[2] = "cs01.s01";
    antNames[3] = "cs01.s02";
    vector<MPosition> antPos(4);
    Vector<double> antDiam(4, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    Vector<double> chanWidth(nchan, 100000);
    Vector<double> chanFreqs(nchan);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
private:
  virtual bool process (const DPBuffer&) { return false; }
  virtual void finish() {}
  virtual void show (std::ostream&) const {}

  int itsNChan, itsNCorr;
};


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
  }
}

void TestPSet::testNone()
{
  TestInput* in = new TestInput(16, 8, 4);
  DPStep::ShPtr step1(in);
  cout << "testNone" << endl;
  ParameterSet parset;
  PreFlagger::PSet pset (in, parset, "");
  pset.updateInfo (in->getInfo());
  ASSERT (!(pset.itsFlagOnBL   || pset.itsFlagOnAmpl || pset.itsFlagOnPhase ||
            pset.itsFlagOnReal || pset.itsFlagOnImag ||
            pset.itsFlagOnAzEl || pset.itsFlagOnUV));
}

void TestPSet::testBL()
{
  TestInput* in = new TestInput(16, 8, 4);
  DPStep::ShPtr step1(in);
  {
    cout << "testBL 1" << endl;
    ParameterSet parset;
    parset.add ("baseline", "[rs01.*, rs02.s01]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (!(pset.itsFlagOnAmpl || pset.itsFlagOnPhase || pset.itsFlagOnReal ||
              pset.itsFlagOnImag || pset.itsFlagOnAzEl  || pset.itsFlagOnUV) &&
            pset.itsFlagOnBL);
    // Make sure the matrix is correct.
    const Matrix<bool>& mat = pset.itsFlagBL;
    ASSERT (mat.shape() == IPosition(2,4,4));
    ASSERT ( mat(0,0) &&  mat(0,1) &&  mat(0,2) &&  mat(0,3));
    ASSERT ( mat(1,0) &&  mat(1,1) &&  mat(1,2) &&  mat(1,3));
    ASSERT ( mat(2,0) &&  mat(2,1) && !mat(2,2) && !mat(2,3));
    ASSERT ( mat(3,0) &&  mat(3,1) && !mat(3,2) && !mat(3,3));
  }
  {
    cout << "testBL 2" << endl;
    ParameterSet parset;
    parset.add ("corrtype", "auto");
    parset.add ("baseline", "[rs01.*, [*s*.*2], rs02.s01]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    // Make sure the matrix is correct.
    const Matrix<bool>& mat = pset.itsFlagBL;
    ASSERT (mat.shape() == IPosition(2,4,4));
    ASSERT ( mat(0,0) && !mat(0,1) && !mat(0,2) && !mat(0,3));
    ASSERT (!mat(1,0) &&  mat(1,1) && !mat(1,2) && !mat(1,3));
    ASSERT (!mat(2,0) && !mat(2,1) && !mat(2,2) && !mat(2,3));
    ASSERT (!mat(3,0) && !mat(3,1) && !mat(3,2) &&  mat(3,3));
  }
  {
    cout << "testBL 3" << endl;
    ParameterSet parset;
    parset.add ("corrtype", "CROSS");
    parset.add ("baseline", "[[rs*, *s*.*1], [cs01.s01,cs01.s02]]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    // Make sure the matrix is correct.
    const Matrix<bool>& mat = pset.itsFlagBL;
    ASSERT (mat.shape() == IPosition(2,4,4));
    ASSERT (!mat(0,0) &&  mat(0,1) &&  mat(0,2) && !mat(0,3));
    ASSERT ( mat(1,0) && !mat(1,1) &&  mat(1,2) && !mat(1,3));
    ASSERT ( mat(2,0) &&  mat(2,1) && !mat(2,2) &&  mat(2,3));
    ASSERT (!mat(3,0) && !mat(3,1) &&  mat(3,2) && !mat(3,3));
  }
  // Some erronous ones.
  cout << "testBL expected error 1" << endl;
  bool err = false;
  try {
    ParameterSet parset;
    parset.add ("corrtype", "crossx");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
  cout << "testBL expected error 2" << endl;
  err = false;
  try {
    ParameterSet parset;
    parset.add ("baseline", "[[a,b,c]]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
  cout << "testBL expected error 3" << endl;
  err = false;
  try {
    ParameterSet parset;
    parset.add ("baseline", "[[a,b], [ ] ]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
}

void TestPSet::testChan()
{
  TestInput* in = new TestInput(16, 32, 4);
  DPStep::ShPtr step1(in);
  {
    cout << "testChan 1" << endl;
    ParameterSet parset;
    parset.add ("chan", "[11..13, 4]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsChannels.size() == 4);
    ASSERT (pset.itsChannels[0] == 4);
    ASSERT (pset.itsChannels[1] == 11);
    ASSERT (pset.itsChannels[2] == 12);
    ASSERT (pset.itsChannels[3] == 13);
    ASSERT (pset.itsChanFlags.shape() == IPosition(2,4,32));
    for (unsigned int i=0; i<32; ++i) {
      if (i==4 || i==11 || i==12 || i==13) {
        ASSERT (allEQ(pset.itsChanFlags.column(i), true));
      } else {
        ASSERT (allEQ(pset.itsChanFlags.column(i), false));
      }
    }
  }
  {
    cout << "testChan 2" << endl;
    ParameterSet parset;
    parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsChannels.size() == 3);
    ASSERT (pset.itsChannels[0] == 1);
    ASSERT (pset.itsChannels[1] == 4);
    ASSERT (pset.itsChannels[2] == 5);
  }
  {
    cout << "testChan 3" << endl;
    ParameterSet parset;
    parset.add ("chan", "[11..13, 4]");
    parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsChannels.size() == 1);
    ASSERT (pset.itsChannels[0] == 4);
  }
}

void TestPSet::testTime()
{
  TestInput* in = new TestInput(16, 8, 4);
  DPStep::ShPtr step1(in);
  {
    cout << "testTime 1" << endl;
    ParameterSet parset;
    parset.add ("abstime", "[1mar2009/12:00:00..2mar2009/13:00:00]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsATimes.size() == 2);
    Quantity q;
    MVTime::read (q, "1mar2009/12:00:00");
    ASSERT (q.getValue("s") == pset.itsATimes[0]);
    ASSERT (pset.itsATimes[1] - pset.itsATimes[0] == 86400+3600);
  }
  {
    cout << "testTime 2" << endl;
    ParameterSet parset;
    parset.add ("reltime", "[12:00:00..13:00:00, 16:00 +- 2min ]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsRTimes.size() == 4);
    ASSERT (pset.itsRTimes[0] == 12*3600);
    ASSERT (pset.itsRTimes[1] == 13*3600);
    ASSERT (pset.itsRTimes[2] == 16*3600-120);
    ASSERT (pset.itsRTimes[3] == 16*3600+120);
  }
  {
    cout << "testTime 3" << endl;
    ParameterSet parset;
    parset.add ("timeofday", "[22:00:00..2:00:00, 23:30 +- 1h ]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsTimes.size() == 8);
    ASSERT (pset.itsTimes[0] == -1);
    ASSERT (pset.itsTimes[1] == 2*3600);
    ASSERT (pset.itsTimes[2] == 22*3600);
    ASSERT (pset.itsTimes[3] == 24*3600+1);
    ASSERT (pset.itsTimes[4] == -1);
    ASSERT (pset.itsTimes[5] == 1800);
    ASSERT (pset.itsTimes[6] == 22*3600+1800);
    ASSERT (pset.itsTimes[7] == 24*3600+1);
  }
  {
    cout << "testTime 4" << endl;
    ParameterSet parset;
    parset.add ("timeslot", "[2..4, 10]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsTimeSlot.size() == 4);
    ASSERT (pset.itsTimeSlot[0] == 2);
    ASSERT (pset.itsTimeSlot[1] == 3);
    ASSERT (pset.itsTimeSlot[2] == 4);
    ASSERT (pset.itsTimeSlot[3] == 10);
  }
  // Some erronous ones.
  cout << "testTime expected error 1" << endl;
  bool err = false;
  try {
    ParameterSet parset;
    parset.add ("reltime", "[12:00:00]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
  cout << "testTime expected error 2" << endl;
  err = false;
  try {
    ParameterSet parset;
    parset.add ("reltime", "[12:00:00..11:00:00]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
  cout << "testTime expected error 3" << endl;
  err = false;
  try {
    ParameterSet parset;
    parset.add ("abstime", "[12:00:00..13:00:00]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
  } catch (std::exception& x) {
    err = true;
    cout << "  " << x.what() << endl; 
  }
  ASSERT (err);
}

void TestPSet::testMinMax()
{
  TestInput* in = new TestInput(16, 8, 4);
  DPStep::ShPtr step1(in);
  {
    cout << "testMinMax 1" << endl;
    ParameterSet parset;
    parset.add ("amplmin", "[23,,,45]");
    parset.add ("amplmax", "112.5");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsFlagOnAmpl);
    ASSERT (pset.itsAmplMin.size() == 4);
    ASSERT (pset.itsAmplMax.size() == 4);
    ASSERT (near(pset.itsAmplMin[0], 23.));
    ASSERT (near(pset.itsAmplMin[1], -1e30));
    ASSERT (near(pset.itsAmplMin[2], -1e30));
    ASSERT (near(pset.itsAmplMin[3], 45.));
    ASSERT (near(pset.itsAmplMax[0], 112.5));
    ASSERT (near(pset.itsAmplMax[1], 112.5));
    ASSERT (near(pset.itsAmplMax[2], 112.5));
    ASSERT (near(pset.itsAmplMax[3], 112.5));
  }
  {
    cout << "testMinMax 2" << endl;
    ParameterSet parset;
    parset.add ("phasemin", "[23]");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsFlagOnPhase);
    ASSERT (pset.itsAmplMin.size() == 4);
    ASSERT (pset.itsAmplMax.size() == 4);
    ASSERT (near(pset.itsPhaseMin[0], 23.));
    ASSERT (near(pset.itsPhaseMin[1], -1e30));
    ASSERT (near(pset.itsPhaseMin[2], -1e30));
    ASSERT (near(pset.itsPhaseMin[3], -1e30));
    ASSERT (near(pset.itsPhaseMax[0], 1e30));
    ASSERT (near(pset.itsPhaseMax[1], 1e30));
    ASSERT (near(pset.itsPhaseMax[2], 1e30));
    ASSERT (near(pset.itsPhaseMax[3], 1e30));
  }
  {
    cout << "testMinMax 3" << endl;
    ParameterSet parset;
    parset.add ("uvmmin", "23");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsFlagOnUV);
    ASSERT (near(pset.itsMinUV, 23.*23.));
    ASSERT (near(pset.itsMaxUV, 1e30));
  }
  {
    cout << "testMinMax 4" << endl;
    ParameterSet parset;
    parset.add ("uvmmax", "23");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsFlagOnUV);
    ASSERT (pset.itsMinUV < 0.);
    ASSERT (near(pset.itsMaxUV, 23.*23.));
  }
  {
    cout << "testMinMax 5" << endl;
    ParameterSet parset;
    parset.add ("uvmmin", "23");
    parset.add ("uvmmax", "123");
    PreFlagger::PSet pset (in, parset, "");
    pset.updateInfo (in->getInfo());
    ASSERT (pset.itsFlagOnUV);
    ASSERT (near(pset.itsMinUV, 23.*23.));
    ASSERT (near(pset.itsMaxUV, 123.*123.));
  }
}

int main()
{
  INIT_LOGGER ("tPSet");
  try {
    TestPSet::testNone();
    TestPSet::testBL();
    TestPSet::testChan();
    TestPSet::testTime();
    TestPSet::testMinMax();
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
