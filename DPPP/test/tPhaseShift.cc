//# tPhaseShift.cc: Test program for class PhaseShift
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
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
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
    info().init (ncorr, nchan, ntime, 0., 10., string());
    MDirection phaseCenter(Quantity(45,"deg"), Quantity(30,"deg"),
                           MDirection::J2000);
    info().set (MPosition(), phaseCenter, phaseCenter, phaseCenter);
    // Define the frequencies.
    Vector<double> chanWidth (nchan, 100000.);
    Vector<double> chanFreqs (nchan);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
    // Fill the baseline stations.
    // Determine nr of stations using:  na*(na+1)/2 = nbl
    // If many baselines, divide into groups of 6 to test if
    // PhaseShift disentangles it correctly.
    int nant = int(-0.5 + sqrt(0.25 + 2*nbl));
    if (nant*(nant+1)/2 < nbl) ++nant;
    int grpszant = 3;
    int grpszbl  = grpszant*(grpszant+1)/2;
    if (nbl > grpszbl) {
      nant = grpszant*(nbl+grpszbl-1)/grpszbl;
    } else {
      grpszant = nant;
      grpszbl  = nbl;
    }
    Vector<Int> ant1(nbl);
    Vector<Int> ant2(nbl);
    int st1 = 0;
    int st2 = 0;
    int lastant = grpszant;
    for (int i=0; i<nbl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (i%grpszbl == grpszbl-1) {
        st1 = lastant;
        st2 = lastant;
        lastant += grpszant;
      } else {
        if (++st2 == lastant) {
          st2 = ++st1;
        }
      }
    }
    Vector<String> antNames(nant);
    vector<MPosition> antPos(nant);
    Vector<double> antDiam(nant, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    itsStatUVW.resize (3, nant);
    for (int i=0; i<nant; ++i) {
      itsStatUVW(0,i) = 0.01 + i*0.02;
      itsStatUVW(1,i) = 0.05 + i*0.03;
      itsStatUVW(2,i) = 0.015 + i*0.025;
    }
  }

  void fillUVW (Matrix<double>& uvw, int count)
  {
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = (itsStatUVW(0,getInfo().getAnt2()[i]) + count*0.002 -
                  (itsStatUVW(0,getInfo().getAnt1()[i]) + count*0.002));
      uvw(1,i) = (itsStatUVW(1,getInfo().getAnt2()[i]) + count*0.004 -
                  (itsStatUVW(1,getInfo().getAnt1()[i]) + count*0.004));
      uvw(2,i) = (itsStatUVW(2,getInfo().getAnt2()[i]) + count*0.006 -
                  (itsStatUVW(2,getInfo().getAnt1()[i]) + count*0.006));
      cout <<getInfo().getAnt1()[i]<<' '<<getInfo().getAnt2()[i]<<' '
           <<uvw(0,i)<<' '<<uvw(1,i)<<' '<<uvw(2,i)<<endl;
    }
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
      data.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
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
    Matrix<double> uvw(3,itsNBl);
    fillUVW (uvw, itsCount);
    buf.setUVW (uvw);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void updateInfo (const DPInfo&) {}
  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
  Matrix<double> itsStatUVW;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput: public DPStep
{
public:
  TestOutput(TestInput* input,
             int ntime, int nbl, int nchan, int ncorr, bool flag)
    : itsInput(input),
      itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    Cube<Complex> result(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(result.size()); ++i) {
      result.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
    }
    Matrix<double> uvw(3,itsNBl);
    itsInput->fillUVW (uvw, itsCount);
    // Check the result.
    ASSERT (allNear(real(buf.getData()), real(result), 1e-7));
    ///cout << imag(buf.getData()) << endl<<imag(result);
    ASSERT (allNear(imag(buf.getData()), imag(result), 1e-7));
    ASSERT (allEQ(buf.getFlags(), itsFlag));
    ASSERT (near(buf.getTime(), 2.+5*itsCount));
    ASSERT (allNear(buf.getUVW(), uvw, 1e-7));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn)
  {
    info() = infoIn;
    MVDirection dir = infoIn.phaseCenter().getValue();
    ASSERT (near(dir.getLong("deg").getValue(), 45.));
    ASSERT (near(dir.getLat("deg").getValue(), 30.));
  }

  TestInput* itsInput;
  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput1: public DPStep
{
public:
  TestOutput1(TestInput* input,
              int ntime, int nbl, int nchan, int ncorr, bool flag)
    : itsInput(input),
      itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    Cube<Complex> result(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(result.size()); ++i) {
      result.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
    }
    Matrix<double> uvw(3,itsNBl);
    itsInput->fillUVW (uvw, itsCount);
    // Check the result.
    ASSERT (! allNear(real(buf.getData()), real(result), 1e-5));
    ASSERT (! allEQ(real(buf.getData()), real(result)));
    ///cout << imag(buf.getData()) << endl<<imag(result);
    ASSERT (! allNear(imag(buf.getData()), imag(result), 1e-5));
    ASSERT (! allEQ(imag(buf.getData()), imag(result)));
    ASSERT (allEQ(buf.getFlags(), itsFlag));
    ASSERT (near(buf.getTime(), 2.+5*itsCount));
    ASSERT (! allNear(buf.getUVW(), uvw, 1e-5));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn)
  {
    info() = infoIn;
    MVDirection dir = infoIn.phaseCenter().getValue();
    ASSERT (near(dir.getLong("deg").getValue(), 50.));
    ASSERT (near(dir.getLat("deg").getValue(), 35.));
  }

  TestInput* itsInput;
  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
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

// Test with a shift to the original center.
void test1(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " flag=" << flag << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  // Keep phase center the same to be able to check if data are correct.
  parset.add ("phasecenter", "[45deg, 30deg]");
  DPStep::ShPtr step2(new PhaseShift(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(in, ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test with a shift to another and then to the original phase center.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " flag=" << flag << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  // First shift to another center, then back to original.
  ParameterSet parset;
  parset.add ("phasecenter", "[50deg, 35deg]");
  ParameterSet parset1;
  parset1.add ("phasecenter", "[]");
  DPStep::ShPtr step2(new PhaseShift(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput1(in, ntime, nbl, nchan, ncorr, flag));
  DPStep::ShPtr step4(new PhaseShift(in, parset1, ""));
  DPStep::ShPtr step5(new TestOutput(in, ntime, nbl, nchan, ncorr, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  step3->setNextStep (step4);
  step4->setNextStep (step5);
  execute (step1);
}


int main()
{
  try {
    test1(10, 3, 32, 4, false);
    test1(10, 10, 30, 1, true);
    test2(10, 6, 32, 4, false);
    test2(10, 6, 30, 1, true);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
