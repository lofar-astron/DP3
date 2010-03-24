//# tPreFlagger.cc: Test program for class PreFlagger
//# Copyright (C) 2010
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/PreFlagger.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>
#include <iostream>

using namespace LOFAR;
using namespace LOFAR::DPPP;
using namespace casa;
using namespace std;

// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set AverageInfo.
  AverageInfo avgInfo;
  DPStep::ShPtr step = step1;
  while (step) {
    step->updateAverageInfo (avgInfo);
    step->show (cout);
    step = step->getNextStep();
  }
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    itsAnt1.resize (nbl);
    itsAnt2.resize (nbl);
    int st1 = 0;
    int st2 = 0;
    for (int i=0; i<nbl; ++i) {
      itsAnt1[i] = st1;
      itsAnt2[i] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    itsAntNames.resize(4);
    itsAntNames[0] = "rs01.s01";
    itsAntNames[1] = "rs02.s01";
    itsAntNames[2] = "cs01.s01";
    itsAntNames[3] = "cs01.s02";
    // Define their positions (more or less WSRT RT0-3).
    itsAntPos.resize (4);
    Vector<double> vals(3);
    vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
    itsAntPos[0] = MPosition(Quantum<Vector<double> >(vals,"m"),
                             MPosition::ITRF);
    vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
    itsAntPos[1] = MPosition(Quantum<Vector<double> >(vals,"m"),
                             MPosition::ITRF);
    vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
    itsAntPos[2] = MPosition(Quantum<Vector<double> >(vals,"m"),
                             MPosition::ITRF);
    vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
    itsAntPos[3] = MPosition(Quantum<Vector<double> >(vals,"m"),
                             MPosition::ITRF);
    // Define the frequencies.
    itsChanFreqs.resize (nchan);
    indgen (itsChanFreqs, 1050000., 100000.);
    // Define start time.
    itsStartTime = 0.5;     // 3 - 0.5*5
  }
private:
  virtual bool process (const DPBuffer&)
  {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-10+itsCount*6);
    }
    DPBuffer buf;
    buf.setTime (itsCount*5 + 3);   //same interval as in updateAveragInfo
    buf.setData (data);
    Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights (weights);
    Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags (flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags (fullResFlags);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
    // Use startchan=0 and timeInterval=5
    { avgInfo.init (itsNCorr, 0, itsNChan, itsNBl, itsNTime, 5); }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of flagged, unaveraged TestInput run by test1/3.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr,
             bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr),
      itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // All baselines except 2-2 should be flagged.
    // Of them only channel 1,4,5 are flagged.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = itsFlag;
    for (int i=0; i<itsNBl; ++i) {
      if (i%16 != 10) {
        for (int j=0; j<itsNChan; ++j) {
          if (j==1 || j==4 || j==5) {
            for (int k=0; k<itsNCorr; ++k) {
              result(k,j,i) = true;
            }
          }
        }
      }
    }
    ///cout << buf.getFlags() << endl << result << endl;
    ASSERT (allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==0);
    ASSERT (int(avgInfo.origNChan())==itsNChan);
    ASSERT (int(avgInfo.nchan())==itsNChan);
    ASSERT (int(avgInfo.ntime())==itsNTime);
    ASSERT (avgInfo.timeInterval()==5);
    ASSERT (int(avgInfo.nchanAvg())==1);
    ASSERT (int(avgInfo.ntimeAvg())==1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
  bool itsFlag;
};

// Test flagging a few antennae and freqs.
void test1(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " flag=" << flag << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add ("baseline", "[rs01.*, *s*.*2, rs02.s01]");
  parset.add ("countflag", "true");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
  step2->showCounts (cout);
}


// Class to check result of flagged, unaveraged TestInput run by test2.
class TestOutput2: public DPStep
{
public:
  TestOutput2(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // A few baselines should be flagged (0, 7, 13, 15)
    // Furthermore channel 1,4,5,11,12,13 are flagged.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = false;
    for (int i=0; i<itsNBl; ++i) {
      if (i%16 == 0  ||  i%16 == 7  ||  i%16 == 13  ||  i%16 == 15) {
        for (int j=0; j<itsNChan; ++j) {
          if (j==4) {
            for (int k=0; k<itsNCorr; ++k) {
              result(k,j,i) = true;
            }
          }
        }
      }
    }
    ///cout << buf.getFlags() << endl << result << endl;
    ASSERT (allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==0);
    ASSERT (int(avgInfo.origNChan())==itsNChan);
    ASSERT (int(avgInfo.nchan())==itsNChan);
    ASSERT (int(avgInfo.ntime())==itsNTime);
    ASSERT (avgInfo.timeInterval()==5);
    ASSERT (int(avgInfo.nchanAvg())==1);
    ASSERT (int(avgInfo.ntimeAvg())==1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
};

// Test flagging a few baselines, freqs, and channels.
void test2(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, false);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add ("chan", "[11..13, 4, 11]");
  parset.add ("baseline", "[[rs01.*,rs01.*],[*s*.*2,*s*.*2],[*s*.*2,rs02.*]]");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test flagging a few antennae or freqs by using multiple steps.
void test3(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test3: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " flag=" << flag << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("expr", "s1");
  parset.add ("s1.freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add ("s1.expr", "s2");
  parset.add ("s1.s2.baseline", "[rs01.*, *s*.*2, rs02.s01]");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}


// Class to check result of flagged, unaveraged TestInput run by test4.
class TestOutput4: public DPStep
{
public:
  TestOutput4(int ntime, int nbl, int nchan, int ncorr,
              bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr),
      itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // All baselines except autocorr should be flagged.
    // Furthermore channel 1,4,5 are flagged.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = true;
    for (int i=0; i<itsNBl; ++i) {
      if (i%16==0 || i%16==5 || i%16==10 || i%16==15) {
        for (int j=0; j<itsNChan; ++j) {
          if (j!=1 && j!=4 && j!=5) {
            for (int k=0; k<itsNCorr; ++k) {
              result(k,j,i) = itsFlag;
            }
          }
        }
      }
    }
    ///cout << buf.getFlags() << endl << result << endl;
    ASSERT (allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==0);
    ASSERT (int(avgInfo.origNChan())==itsNChan);
    ASSERT (int(avgInfo.nchan())==itsNChan);
    ASSERT (int(avgInfo.ntime())==itsNTime);
    ASSERT (avgInfo.timeInterval()==5);
    ASSERT (int(avgInfo.nchanAvg())==1);
    ASSERT (int(avgInfo.ntimeAvg())==1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
  bool itsFlag;
};

// Test flagging a few antennae and freqs by using multiple steps.
void test4(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test4: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " flag=" << flag << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("expr", "(s1&s1),(s2|s2)");
  parset.add ("s1.freqrange", "[ 1.1 .. 1.2 MHz, 1.5MHz+-65000Hz]");
  parset.add ("s2.baseline", "[rs01.*, *s*.*2, rs02.s01]");
  parset.add ("s2.corrtype", "cross");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput4(ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

typedef bool CheckFunc (Complex value, double time, int ant1, int ant2);

// Class to check result of flagged, unaveraged TestInput run by test5.
class TestOutput5: public DPStep
{
public:
  TestOutput5 (CheckFunc* cfunc)
    : itsCount(0),
      itsCFunc (cfunc)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    const Cube<Complex>& data = buf.getData();
    const IPosition& shp = data.shape();
    Cube<bool> result(shp);
    for (int i=0; i<shp[2]; ++i) {
      int a1 = i/4;
      int a2 = i%4;
      for (int j=0; j<shp[1]; ++j) {
        bool flag = false;
        for (int k=0; k<shp[0]; ++k) {
          if (!flag) flag = itsCFunc(data(k,j,i), buf.getTime(), a1, a2);
        }
        for (int k=0; k<shp[0]; ++k) {
          result(k,j,i) = flag;
        }
      }
    }
    ///cout << buf.getFlags() << endl << result << endl;
    ASSERT (allEQ(buf.getFlags(), result));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo&) {}

  int itsCount;
  CheckFunc* itsCFunc;
};

// Test flagging on a single parameter.
void test5(const string& key, const string& value, CheckFunc* cfunc)
{
  cout << "test5: " << key << '=' << value << endl;
  // Create the steps.
  TestInput* in = new TestInput(2, 6, 5, 4, false);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add (key, value);
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput5(cfunc));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test flagging on multiple parameters.
void test6(const string& key1, const string& value1,
           const string& key2, const string& value2, CheckFunc* cfunc)
{
  cout << "test6: " << key1 << '=' << value1 << ' '
       << key2 << '=' << value2 << endl;
  // Create the steps.
  TestInput* in = new TestInput(6, 10, 8, 4, false);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add (key1, value1);
  parset.add (key2, value2);
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput5(cfunc));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

bool checkBL(Complex, double, int a1, int a2)
  { return a1==a2; }
bool checkAmplMin (Complex data, double, int, int)
  { return abs(data) < 9.5; }
bool checkAmplMax (Complex data, double, int, int)
  { return abs(data) > 31.5; }
bool checkPhaseMin (Complex data, double, int, int)
  { return arg(data) < 1.4; }
bool checkPhaseMax (Complex data, double, int, int)
  { return arg(data) > 2.1; }
bool checkRealMin (Complex data, double, int, int)
  { return real(data) < 5.5; }
bool checkRealMax (Complex data, double, int, int)
  { return real(data) > 29.4; }
bool checkImagMin (Complex data, double, int, int)
  { return imag(data) < -1.4; }
bool checkImagMax (Complex data, double, int, int)
  { return imag(data) > 20.5; }
bool checkTimeSlot (Complex, double time, int, int)
  { return time<5; }
bool checkNone (Complex, double, int, int)
  { return false; }
bool checkAll (Complex, double, int, int)
  { return true; }

// Test flagging on various fields.
void testMany()
{
  test5("corrtype", "auto", &checkBL);
  test5("amplmin", "9.5", &checkAmplMin);
  test5("amplmax", "31.5", &checkAmplMax);
  test5("phasemin", "1.4", &checkPhaseMin);
  test5("phasemax", "2.1", &checkPhaseMax);
  test5("realmin", "5.5", &checkRealMin);
  test5("realmax", "29.4", &checkRealMax);
  test5("imagmin", "-1.4", &checkImagMin);
  test5("imagmax", "20.5", &checkImagMax);
  test5("timeslot", "0", &checkTimeSlot);
  test5("abstime", "17-nov-1858/0:0:2..17nov1858/0:0:4", &checkTimeSlot);
  test5("abstime", "17-nov-1858/0:0:3+-1s", &checkTimeSlot);
  test5("reltime", "0:0:0..0:0:6", &checkTimeSlot);
  test5("reltime", "0:0:2+-1s", &checkTimeSlot);
  test5("reltime", "0:0:2+-20s", &checkAll);
  test6("abstime", "17-nov-1858/0:0:20+-19s",
        "reltime", "0:0:20+-20s", &checkAll);
  test6("abstime", "17-nov-1858/0:0:3+-2s",
        "reltime", "0:0:20+-20s", &checkTimeSlot);
  test6("abstime", "17-nov-1858/0:0:20+-19s",
        "reltime", "0:0:3+-2s", &checkTimeSlot);
  test6("abstime", "17-nov-1858/0:0:20+-9s",
        "reltime", "0:0:3+-2s", &checkNone);
  // Elevation is 12738s; azimuth=86121s
  test5("elevation", "180deg..190deg", &checkNone);
  test5("elevation", "12730s..12740s", &checkAll);
  test5("azimuth", "180deg..190deg", &checkNone);
  test5("azimuth", "86120s..86125s", &checkAll);
  test6("azimuth", "86120s..86125s", "elevation", "180deg..190deg", &checkNone);
  test6("azimuth", "86120s..86125s", "elevation", "12730s..12740s", &checkAll);
  test5("lst", "0.154d..0.155d", &checkAll);
}

int main()
{
  INIT_LOGGER ("tPreFlagger");
  try {

    test1(10, 16, 32, 4, false);
    test1(10, 16, 32, 4, true);
    test2( 2, 16, 32, 4);
    test2( 2, 36, 16, 2);
    test3( 3, 16, 32, 4, false);
    test3( 4, 16,  4, 2, true);
    test4( 3, 16, 32, 4, false);
    testMany();
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
