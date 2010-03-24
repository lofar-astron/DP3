//# tUVWFlagger.cc: Test program for class UVWFlagger
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
#include <DPPP/UVWFlagger.h>
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
  TestInput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
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
    indgen (itsChanFreqs, 10500000., 1000000.);
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
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime (itsCount*30 + 4472025740.0);
    buf.setData (data);
    buf.setUVW  (uvw);
    Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags (flags);
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
};

// Class to check result of flagged, unaveraged TestInput run by test1.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // Flag where u,v,w matches intervals given in test1.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = false;
    for (int i=0; i<itsNBl; ++i) {
      double u = 1+i+itsCount;
      double v = 2+i+itsCount;
      double w = 3+i+itsCount;
      double uv = sqrt(u*u+v*v);
      if ((uv>5.5 && uv<8.5) || (u>20.5 && u<23.5) || (u>31.5 && u<40.5)
          || (v>11.5 && v<14.5) || w<3.5 || w>44.5) {
        for (int j=0; j<itsNChan; ++j) {
          for (int k=0; k<itsNCorr; ++k) {
            result(k,j,i) = true;
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
    // Flag where u,v,w matches intervals given in test1.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = false;
    for (int i=0; i<itsNBl; ++i) {
      for (int j=0; j<itsNChan; ++j) {
        double wavel = 2.99792458e+08 / (10.5e6 + j*1e6);
        double u = (1+i+itsCount) / wavel;
        double v = (2+i+itsCount) / wavel;
        double w = (3+i+itsCount) / wavel;
        double uv = sqrt(u*u+v*v);
        if ((uv>0.2 && uv<0.31) || (u>1.55 && u<1.485) || (u>0.752 && u<0.862)
            || (v>0.42 && v<0.53) || w<0.12 || w>1.63) {
          for (int k=0; k<itsNCorr; ++k) {
            result(k,j,i) = true;
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


// Class to check result of flagged, unaveraged TestInput run by test3.
class TestOutput3: public DPStep
{
public:
  TestOutput3(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {
    ASSERT (ntime==2 && nbl==16);
  }
private:
  virtual bool process (const DPBuffer& buf)
  {
    // These are the UVW coordinates as calculated by UVWFlagger for the
    // station positions and phase center defined in TestInput.
    double uvwvals[] = {
      0, 0, 0,
      0.423756, -127.372, 67.1947,
      0.847513, -254.744, 134.389,
      0.277918, -382.015, 201.531,
      -0.423756, 127.372, -67.1947,
      0, 0, 0,
      0.423756, -127.372, 67.1947,
      -0.145838, -254.642, 134.336,
      -0.847513, 254.744, -134.389,
      -0.423756, 127.372, -67.1947,
      0, 0, 0,
      -0.569594, -127.27, 67.1417,
      -0.277918, 382.015, -201.531,
      0.145838, 254.642, -134.336,
      0.569594, 127.27, -67.1417,
      0, 0, 0,
      0, 0, 0,
      0.738788, -127.371, 67.1942,
      1.47758, -254.742, 134.388,
      1.22276, -382.013, 201.53,
      -0.738788, 127.371, -67.1942,
      0, 0, 0,
      0.738788, -127.371, 67.1942,
      0.483976, -254.642, 134.336,
      -1.47758, 254.742, -134.388,
      -0.738788, 127.371, -67.1942,
      0, 0, 0,
      -0.254812, -127.271, 67.1421,
      -1.22276, 382.013, -201.53,
      -0.483976, 254.642, -134.336,
      0.254812, 127.271, -67.1421,
      0, 0, 0
    };
    Cube<double> uvws(IPosition(3,3,16,2), uvwvals, SHARE);
    // Flag where u,v,w matches intervals given in test3.
    Cube<bool> result(itsNCorr,itsNChan,itsNBl);
    result = false;
    for (int i=0; i<itsNBl; ++i) {
      double u = uvws(0,i,itsCount);
      double v = uvws(1,i,itsCount);
      double w = uvws(2,i,itsCount);
      double uv = sqrt(u*u+v*v);
      if ((uv>5.5 && uv<8.5) || (u>20.5 && u<23.5) || (u>31.5 && u<40.5)
          || (v>11.5 && v<14.5) || w<3.5 || w>44.5) {
        for (int j=0; j<itsNChan; ++j) {
          for (int k=0; k<itsNCorr; ++k) {
            result(k,j,i) = true;
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

// Test flagging a few baselines on UV in m.
void test1(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("uvmrange", "[5.5..8.5]");
  parset.add ("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add ("vmrange", "[11.5..14.5]");
  parset.add ("wmmax", "44.5");
  parset.add ("wmmin", "3.5");
  DPStep::ShPtr step2(new UVWFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test flagging a few baselines on UV in wavelengths.
void test2(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("uvlambdarange", "[0.2..0.31]");
  parset.add ("ulambdarange", "[1.55..1.485, 0.807+-0.055]");
  parset.add ("vlambdarange", "[0.42..0.53]");
  parset.add ("wlambdamax", "1.63");
  parset.add ("wlambdamin", "0.12");
  DPStep::ShPtr step2(new UVWFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
  step2->showCounts (cout);
}

// Test flagging a few baselines on UV in m with a different phase center.
void test3(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test3: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("uvmrange", "[5.5..8.5]");
  parset.add ("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add ("vmrange", "[11.5..14.5]");
  parset.add ("wmmax", "44.5");
  parset.add ("wmmin", "3.5");
  parset.add ("phasecenter", "[-1.92653768rad, 1.09220917rad, j2000]");
  DPStep::ShPtr step2(new UVWFlagger(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput3(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test constructing with the Sun as phase center.
void test4()
{
  cout << "test4" << endl;
  // Create the steps.
  TestInput* in = new TestInput(1,1,1,1);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("uvmrange", "[5.5..8.5]");
  parset.add ("phasecenter", "Sun");
  DPStep::ShPtr step2(new UVWFlagger(in, parset, ""));
  step2->show (cout);
}


int main()
{
  INIT_LOGGER ("tUVWFlagger");
  try {

    test1( 10,  16, 32, 4);
    test1(100, 105, 32, 4);
    test2(  2,  16, 32, 4);
    test2(  2,  36, 16, 2);
    test2( 10,  16, 32, 4);
    test2(100, 105, 32, 4);
    test3(  2,  16, 32, 4);
    test4();
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
