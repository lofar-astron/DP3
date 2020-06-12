// tInterpolateStep.cc: Test program for class InterpolateStep
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
// $Id: tInterpolateStep.cc 24221 2013-03-12 12:24:48Z diepen $
//
// @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP_Interpolate/Interpolate.h>
#include <DPPP/DPRun.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>
#include <iostream>

using namespace LOFAR;
using namespace DP3::DPPP;
using namespace casa;
using namespace std;

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nant, int nchan, int ncorr, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nant*(nant+1)/2), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    Vector<Int> ant1(itsNBl);
    Vector<Int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (int i=0; i<itsNBl; ++i) {
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
    vector<MPosition> antPos (4);
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
    Vector<double> chanFreqs(nchan);
    Vector<double> chanWidth(nchan, 100000.);
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
    cout << "Input step " << itsCount << ' '<< itsCount*5+2<<endl;
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(1.6, 0.9);
    }
    if (itsCount == 5) {
      data += Complex(10.,10.);
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
  virtual void updateInfo (const DPInfo&)
    // Use startchan=0 and timeInterval=5
    { info().init (itsNCorr, itsNChan, itsNTime, 100, 5, string(), string()); }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nant, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nant*(nant+1)/2), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    cout << "Output step " << itsCount << ' '<<itsCount*5+2<<endl;
    // Fill expected result in similar way as TestInput.
    Cube<Complex> result(itsNCorr,itsNChan,itsNBl);
    for (int i=0; i<int(result.size()); ++i) {
      result.data()[i] = Complex(1.6, 0.9);
    }
    if (itsCount == 5) {
      result += Complex(10.,10.);
    }
    // Check the result.
    ///cout << buf.getData()<< result;
    ASSERT (allNear(real(buf.getData()), real(result), 1e-10));
    ASSERT (allNear(imag(buf.getData()), imag(result), 1e-10));
    ASSERT (near(buf.getTime(), 2+5.*itsCount));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    ASSERT (int(info.origNChan())==itsNChan);
    ASSERT (int(info.nchan())==itsNChan);
    ASSERT (int(info.ntime())==itsNTime);
    ASSERT (info.startTime()==100);
    ASSERT (info.timeInterval()==5);
    ASSERT (int(info.nchanAvg())==1);
    ASSERT (int(info.ntimeAvg())==1);
    ASSERT (int(info.chanFreqs().size()) == itsNChan);
    ASSERT (int(info.chanWidths().size()) == itsNChan);
    ASSERT (info.msName().empty());
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set the info in each step.
  step1->setInfo (DPInfo());
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
  DPStep::ShPtr step = step1;
  while (step) {
    step->showCounts (cout);
    step = step->getNextStep();
  }
}

void test1(int ntime, int nant, int nchan, int ncorr, bool flag, int threshold)
{
  cout << "test1: ntime=" << ntime << " nrant=" << nant << " nchan=" << nchan
       << " ncorr=" << ncorr << " threshold=" << threshold << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("windowsize", "9");
  DPStep::ShPtr step2 = DPRun::findStepCtor("interpolate")(in, parset, "");
  DPStep::ShPtr step3(new TestOutput(ntime, nant, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  step2->show (cout);
  execute (step1);
}

int main()
{
  try {
    test1(10, 2, 32, 4, false, 1);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
