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
    itsChanFreqs.resize (10);
    indgen (itsChanFreqs, 1050000., 100000.);
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
    buf.setTime (itsCount*5 + 2);   //same interval as in updateAveragInfo
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

// Class to check result of flagged, unaveraged TestInput run by test1.
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
    // All baselines except 22 should be flagged.
    // Furthermore channel 1,4,5 are flagged.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = true;
    for (int i=10; i<itsNBl; i+=16) {
      for (int j=0; j<itsNChan; ++j) {
        if (j!=1 && j!=4 && j!=5) {
          for (int k=0; k<itsNCorr; ++k) {
            result(k,j,i) = itsFlag;
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
    // A few baselines should be flagged (0, 13, 15, ...)
    // Furthermore channel 1,4,5,11,12,13 are flagged.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = true;
    for (int i=0; i<itsNBl; ++i) {
      if (i%16 != 0  &&  i%16 != 13  &&  i%16 != 15) {
        for (int j=0; j<itsNChan; ++j) {
          if (j!=1 && j!=4 && j!=5 && j!=11 && j!=12 && j!=13) {
            for (int k=0; k<itsNCorr; ++k) {
              result(k,j,i) = false;
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
  parset.add ("antenna", "[rs01.*, *s*.*2, rs02.s01]");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

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
  parset.add ("antenna1", "[rs01.*, *s*.*2, *s*.*2]");
  parset.add ("antenna2", "[rs01.*, *s*.*2, rs02.*]");
  DPStep::ShPtr step2(new PreFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}


int main()
{
  INIT_LOGGER ("tPreFlagger");
  try {

    test1(10, 16, 32, 4, false);
    test1(10, 16, 32, 4, true);
    test2( 2, 16, 32, 4);
    test2( 2, 36, 16, 2);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
