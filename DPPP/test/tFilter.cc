//# tFilter.cc: Test program for class Filter
//# Copyright (C) 2012
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
#include <DPPP/Filter.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
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
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {
    // Define start time 0.5 (= 3 - 0.5*5) and time interval 5.
    info().init (ncorr, nchan, ntime, 0.5, 5., string());
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
    // Define their positions (more or less WSRT RT0-3).
    vector<MPosition> antPos(4);
    Vector<double> vals(3);
    vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
    antPos[0] = MPosition(Quantum<Vector<double> >(vals,"m"),
                          MPosition::ITRF);
    vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
    antPos[1] = MPosition(Quantum<Vector<double> >(vals,"m"),
                          MPosition::ITRF);
    vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
    antPos[2] = MPosition(Quantum<Vector<double> >(vals,"m"),
                          MPosition::ITRF);
    vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
    antPos[3] = MPosition(Quantum<Vector<double> >(vals,"m"),
                          MPosition::ITRF);
    Vector<double> antDiam(4, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    Vector<double> chanWidth(nchan, 100000.);
    Vector<double> chanFreqs(nchan);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
private:
  virtual bool process (const DPBuffer&)
  {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    DPBuffer buf;
    buf.setTime (itsCount*5 + 2);
    buf.setExposure (0.1*(itsCount+1));
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
    }
    buf.setData (data);
    Cube<float> weights(data.shape());
    buf.setWeights (weights);
    indgen (weights);
    Cube<bool> flags(data.shape());
    flags = itsFlag;
    // Set part of the flags to another value.
    flags(IPosition(3,0), flags.shape()-1, IPosition(3,1,3,4)) = !itsFlag;
    buf.setFlags (flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // Assume they are averaged for 2 chan, 2 time.
    Cube<bool> fullResFlags;
    if (itsNCorr == 4) {
      fullResFlags = flags.copy().reform(IPosition(3,2*itsNChan,2,itsNBl));
    } else {
      fullResFlags.resize (IPosition(3, itsNChan, 1, itsNBl));
      fullResFlags = true;
    }
    buf.setFullResFlags (fullResFlags);
    Matrix<double> uvw(3,itsNBl);
    indgen (uvw, double(itsCount*100));
    buf.setUVW (uvw);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo&)
  {
    // Use timeInterval=5
    info().init (itsNCorr, itsNChan, itsNTime, 100, 5, string());
    // Define the frequencies.
    Vector<double> chanFreqs(itsNChan);
    Vector<double> chanWidth(itsNChan, 100000.);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of averaging TestInput.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr,
             int nblout, int stchan, int nchanOut, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsNBlOut(nblout),
      itsStChan(stchan), itsNChanOut(nchanOut),
      itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // Fill expected result in similar way as TestInput.
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
    }
    Cube<float> weights(data.shape());
    indgen (weights);
    Cube<bool> flags(data.shape());
    flags = itsFlag;
    // Set part of the flags to another value.
    flags(IPosition(3,0), flags.shape()-1, IPosition(3,1,3,4)) = !itsFlag;
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // Assume they are averaged for 2 chan, 2 time.
    Cube<bool> fullResFlags;
    if (itsNCorr == 4) {
      fullResFlags = flags.copy().reform(IPosition(3,2*itsNChan,2,itsNBl));
    } else {
      fullResFlags.resize (IPosition(3, itsNChan, 1, itsNBl));
      fullResFlags = true;
    }
    Matrix<double> uvw(3,itsNBl);
    indgen (uvw, double(itsCount*100));
    Slicer slicer(IPosition(3,0,itsStChan,0),
                  IPosition(3,itsNCorr,itsNChanOut,itsNBlOut));
    // Check the expected result.
    ASSERT (allEQ(buf.getData(), data(slicer)));
    ASSERT (allEQ(buf.getFlags(), flags(slicer)));
    ASSERT (allEQ(buf.getWeights(), weights(slicer)));
    ASSERT (allEQ(buf.getUVW(), uvw(IPosition(2,0,0),
                                    IPosition(2,2,itsNBlOut-1))));
    if (itsNCorr == 4) {
      ASSERT (allEQ(buf.getFullResFlags(),
                    fullResFlags(Slicer(IPosition(3,itsStChan*2,0,0),
                                        IPosition(3,2*itsNChanOut,2,itsNBlOut)))));
    } else {
      ASSERT (allEQ(buf.getFullResFlags(),
                    fullResFlags(Slicer(IPosition(3,itsStChan,0,0),
                                        IPosition(3,itsNChanOut,1,itsNBlOut)))));
    }
    ASSERT (near(buf.getTime(), itsCount*5.+2));
    ASSERT (near(buf.getExposure(), 0.1*(itsCount+1)));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    ASSERT (int(info.origNChan())==itsNChan);
    ASSERT (int(info.nchan())==itsNChanOut);
    ASSERT (int(info.nbaselines())==itsNBlOut);
    ASSERT (int(info.ntime())==itsNTime);
    ASSERT (info.timeInterval()==5.);
    ASSERT (int(info.nchanAvg())==1);
    ASSERT (int(info.ntimeAvg())==1);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNBlOut, itsStChan, itsNChanOut;
  bool itsFlag;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set DPInfo.
  step1->setInfo (DPInfo());
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test filtering of channels only.
void test1(int ntime, int nbl, int nchan, int ncorr,
           int startchan, int nchanout, bool flag)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " startchan=" << startchan
       << " nchanout=" << nchanout << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("startchan", toString(startchan));
  parset.add ("nchan", toString(nchanout)+"+nchan-nchan");
  DPStep::ShPtr step2(new Filter(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     nbl, startchan, nchanout, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test filtering of baselines and channels.
void test2(int ntime, int nbl, int nchan, int ncorr,
           int startchan, int nchanout, bool flag)
{
  ASSERT (nbl<=4); // otherwise baseline selection removes more than the first
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " startchan=" << startchan
       << " nchanout=" << nchanout << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("startchan", toString(startchan)+"+nchan-nchan");
  parset.add ("nchan", toString(nchanout));
  // This removes the first baseline.
  parset.add ("baseline", "[[rs01.s01,rs*]]");
  DPStep::ShPtr step2(new Filter(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     2, startchan, nchanout, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}


int main()
{
  try {
    test1(10,  3, 32, 4, 2, 24, false);
    test1(10, 10, 30, 1, 3,  3, true);
    test1(10, 10,  1, 4, 0,  1, true);
    test2(10,  4, 32, 4, 2, 24, false);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
