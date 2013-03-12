//# tScaleData.cc: Test program for class ScaleData
//# Copyright (C) 2013
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
//# $Id: tScaleData.cc 23691 2013-02-12 13:32:33Z diepen $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/ScaleData.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/StreamUtil.h>
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
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {
    info().init (ncorr, nchan, ntime, 0., 5., string(), string());
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
    Vector<double> chanWidth(nchan, 1000000.);
    Vector<double> chanFreqs(nchan);
    indgen (chanFreqs, 10500000., 1000000.);
    info().set (chanFreqs, chanWidth);
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
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 0.5f, 0.01f);
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime (itsCount*30 + 4472025740.0);
    buf.setData (data);
    buf.setWeights (weights);
    buf.setUVW  (uvw);
    Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags (flags);
    Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = false;
    buf.setFullResFlags (fullResFlags);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo&) {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of TestInput run by test1.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  void addData (Cube<Complex>& to, const Cube<Complex>& from, int bl)
  {
    to += from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl));
  }
  void addConjData (Cube<Complex>& to, const Cube<Complex>& from, int bl)
  {
    to += conj(from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl)));
  }
  virtual bool process (const DPBuffer& buf)
  {
    // Fill data and scale as needed.
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    Complex* dataPtr = data.data();
    int cnt=0;
    for (int i=0; i<itsNBl; ++i) {
      double freq = 10.5;
      for (int j=0; j<itsNChan; ++j) {
        double sc1 = 3 + 2*freq + freq*freq;
        double sc2 = sc1;
        if ((i%16)/4 == 0) {
          sc1 = 2 + 0.5*freq;
        }
        if (i%4 == 0) {
          sc2 = 2 + 0.5*freq;
        }
        double scale = sqrt(sc1*sc2);
        freq += 1;
        for (int k=0; k<itsNCorr; ++k) {
          *dataPtr++ = Complex(cnt+itsCount*10,cnt-10+itsCount*6) * scale;
          cnt++;
        }
      }
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 0.5f, 0.01f);
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    ASSERT (allEQ (buf.getData(), data));
    ASSERT (buf.getFlags().shape() == IPosition(3,itsNCorr,itsNChan,itsNBl));
    ASSERT (allEQ (buf.getFlags(), false));
    ASSERT (allEQ (buf.getWeights(), weights));
    ASSERT (allEQ (buf.getUVW(), uvw));
    ASSERT (buf.getFullResFlags().shape() == IPosition(3,itsNChan,1,itsNBl));
    ASSERT (allEQ (buf.getFullResFlags(), false));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn)
  {
    info() = infoIn;
    ASSERT (int(infoIn.origNChan())==itsNChan);
    ASSERT (int(infoIn.nchan())==itsNChan);
    ASSERT (int(infoIn.ntime())==itsNTime);
    ASSERT (infoIn.timeInterval()==5);
    ASSERT (int(infoIn.nchanAvg())==1);
    ASSERT (int(infoIn.ntimeAvg())==1);
    ASSERT (int(infoIn.nbaselines())==itsNBl);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set DPInfo.
  step1->setInfo (DPInfo());
  step1->getNextStep()->show (cout);
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test scaling.
void test1(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("stations", "[rs01.s01, *]");
  parset.add ("coeffs", "[[2,0.5],[3,2,1]]");
  parset.add ("scalesize", "false");
  DPStep::ShPtr step2(new ScaleData(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}


int main()
{
  INIT_LOGGER ("tScaleData");
  try {
    test1 ( 2,  4, 4, 1);
    test1 (10,  16, 32, 4);
    test1 (10,  12, 16, 2);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
