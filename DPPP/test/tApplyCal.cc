//# tApplyCal.cc: Test program for class AORFlagger
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
//# $Id: tApplyCal.cc 24221 2013-08-02 12:24:48Z tammo $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/ApplyCal.h>
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
// 9 baselines, 3 antennas, 4 correlations
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nchan)
    : itsCount(0), itsNTime(ntime), itsNChan(nchan), itsNBl(9), itsNCorr(4),
      itsTimeInterval(5.)
  {
    info().init (itsNCorr, nchan, ntime, 4472025740.0, itsTimeInterval,
        string(), string());
    // Fill the baseline stations; use 3 stations.
    // So they are called 00 01 02 10 11 12 20 21 22, etc.

    Vector<Int> ant1(itsNBl);
    Vector<Int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (int i=0; i<itsNBl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 3) {
        st2 = 0;
        if (++st1 == 3) {
          st1 = 0;
        }
      }
    }
    Vector<String> antNames(4);
    antNames[0] = "ant1";
    antNames[1] = "ant2";
    antNames[2] = "ant3";
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
      data.data()[i] = Complex(1,0);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    weights=1.;

    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime (itsCount*itsTimeInterval + 4472025740.0);
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

  int itsCount, itsNTime, itsNChan, itsNBl, itsNCorr, itsTimeInterval;
};



// Class to check result of TestInput run by tests.
class TestOutput: public DPStep
{
public:
  enum tests {WeightsNotChanged=1, DataNotChanged=2, DataChanged=4,
  DataEquals=8, WeightEquals=16};
  TestOutput(int ntime, int nchan, int doTest)
    : itsCount(0), itsTimeStep(0), itsNTime(ntime), itsNBl(9), itsNChan(nchan),
      itsNCorr(4), itsTimeInterval(5.), itsDoTest(doTest)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    // Fill data and scale as needed.
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);

    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(1,0);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 1.0f, 0.0f);

    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }

    // The same gain corrections as in tApplyCal_tmp.parmdb
    vector<Cube<Complex> > gains(4); // cube for every corr
    for (int corr=0; corr<4; ++corr) {
      gains[corr].resize(IPosition(3,2,2,3)); //freq,time,ant;
    }

    gains[0](Slice(0,2),Slice(0,2),Slice(0,3))=1;
    gains[1](Slice(0,2),Slice(0,2),Slice(0,3))=0;
    gains[2](Slice(0,2),Slice(0,2),Slice(0,3))=0;
    gains[3](Slice(0,2),Slice(0,2),Slice(0,3))=1;
    // ant2
    gains[0](0,0,1)=2;
    gains[3](0,0,1)=3;
    gains[0](1,1,1)=Complex(3.,4.);
    // ant3
    gains[2](1,0,2)=.5;

    if (itsDoTest & DataEquals) {
      //cout<<"weights="<<buf.getWeights()<<endl;
      for (int bl=0; bl<itsNBl; ++bl) {
        for (int chan=0; chan<itsNChan; ++chan) {
          for (int corr=0; corr<itsNCorr; ++corr) {
            data(corr,chan,bl) /=
              (gains[corr/2*3](chan/(itsNChan/2),
                  itsTimeStep/(itsNTime/2), info().getAnt1()[bl]) *
               conj(gains[corr%2*3](chan/(itsNChan/2),
                  itsTimeStep/(itsNTime/2), info().getAnt2()[bl])));
          }
          if (info().getAnt2()[bl]==2 && (itsTimeStep/(itsNTime/2))==0
              && (chan/(itsNChan/2))==1) {
            data(0,chan,bl)-=0.5*data(1,chan,bl);
            data(2,chan,bl)-=0.5*data(3,chan,bl);
          }
          if (info().getAnt1()[bl]==2 && (itsTimeStep/(itsNTime/2))==0
              && (chan/(itsNChan/2))==1) {
            data(0,chan,bl)-=0.5*data(2,chan,bl);
            data(1,chan,bl)-=0.5*data(3,chan,bl);
          }
        }
      }

    }

    if (itsDoTest & WeightEquals) {
      ASSERT ( near(buf.getWeights()(0,0,1),4.));
      ASSERT ( near(buf.getWeights()(1,0,1),9.));
      ASSERT ( near(buf.getWeights()(2,0,1),4.));
      ASSERT ( near(buf.getWeights()(3,0,1),9.));
      ASSERT ( near(buf.getWeights()(0,31,5),0.8));
    }

    if (itsDoTest & DataEquals) {
      ASSERT (allNear (buf.getData(), data, 1.e-7));
    }

    if (itsDoTest & DataNotChanged) {
      ASSERT (allNear (buf.getData(), data, 1.e-7));
    }
    if (itsDoTest & DataChanged) {
      ASSERT (!(allNear (buf.getData(), data, 1.e-7)));
    }
    if (itsDoTest & WeightsNotChanged) {
      ASSERT (allNear (buf.getWeights(), weights, 1.e-7));
    }
    itsCount++;
    itsTimeStep++;
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
    ASSERT (infoIn.timeInterval()==itsTimeInterval);
    ASSERT (int(infoIn.nbaselines())==itsNBl);
  }

  int itsCount;
  int itsTimeStep;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsTimeInterval, itsDoTest;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set DPInfo.
  step1->setInfo (DPInfo());

  const DPStep::ShPtr& step=step1->getNextStep();

  // TODO: do line below for any step that is an ApplyCal
  step->show (cout);

  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test clock + tec, and test two ApplyCals in sequence
void testclocktec(int ntime, int nchan)
{
  cout << "testclocktec: ntime=" << ntime << " nchan=" << nchan << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nchan);
  DPStep::ShPtr step1(in);

  ParameterSet parset1;
  parset1.add ("correction", "tec");
  parset1.add ("parmdb", "tApplyCal_tmp.parmdb");
  parset1.add ("timeslotsperparmupdate", "5");
  parset1.add ("updateweights", "true");
  DPStep::ShPtr step2(new ApplyCal(in, parset1, ""));

  ParameterSet parset2;
  parset2.add ("correction", "clock");
  parset2.add ("parmdb", "tApplyCal_tmp.parmdb");
  parset2.add ("timeslotsperparmupdate", "5");
  parset2.add ("updateweights", "true");
  DPStep::ShPtr step3(new ApplyCal(in, parset2, ""));

  ParameterSet parset3;
  parset3.add ("correction", "commonscalarphase");
  parset3.add ("parmdb", "tApplyCal_tmp.parmdb");
  parset3.add ("timeslotsperparmupdate", "1");
  parset3.add ("udpateweights", "true");
  DPStep::ShPtr step4(new ApplyCal(in, parset3, ""));

  DPStep::ShPtr step5(new TestOutput(ntime, nchan,
      TestOutput::DataChanged | TestOutput::WeightsNotChanged));

  step1->setNextStep (step2);
  step2->setNextStep (step3);
  step3->setNextStep (step4);
  step4->setNextStep (step5);
  execute (step1);
  cout<<"finished testclocktec"<<endl;
}

// Test gain
void testgain(int ntime, int nchan)
{
  cout << "testgain: ntime=" << ntime << " nchan=" << nchan << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nchan);
  DPStep::ShPtr step1(in);

  ParameterSet parset1;
  parset1.add ("correction", "gain");
  parset1.add ("parmdb", "tApplyCal_tmp.parmdb");
  parset1.add ("timeslotsperparmupdate", "5");
  parset1.add ("updateweights", "true");
  DPStep::ShPtr step2(new ApplyCal(in, parset1, ""));

  DPStep::ShPtr step3(new TestOutput(ntime, nchan,
      TestOutput::DataEquals));

  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}


int main()
{
  INIT_LOGGER ("tApplyCal");
  try {
    testclocktec (10,  32);
    testgain (10, 32);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
   return 1;
  }
  return 0;
}
