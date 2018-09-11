//# tApplyCalH5.cc: Test program for class ApplyCal
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
//# $Id: tApplyCalH5.cc 24221 2013-08-02 12:24:48Z tammo $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/H5Parm.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/StreamUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>
#include <iostream>

using namespace LOFAR;
using namespace DP3::DPPP;
using namespace casa;
using namespace std;

// Simple class to generate input arrays.
// 9 baselines, 3 antennas, 4 correlations
class TestInput: public DPInput
{
public:
  TestInput(uint ntime, uint nchan)
    : itsCount(0), itsNTime(ntime), itsNChan(nchan), itsNBl(9), itsNCorr(4),
      itsTimeInterval(5.), itsFirstTime(4472025740.0)
  {
    info().init (itsNCorr, nchan, ntime, itsFirstTime, itsTimeInterval,
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
    Vector<String> antNames(3);
    antNames[0] = "ant1";
    antNames[1] = "ant2";
    antNames[2] = "ant3";
    // Define their positions (more or less WSRT RT0-3).
    vector<MPosition> antPos(3);
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
    Vector<double> antDiam(3, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    Vector<double> chanWidth(nchan, 100.e6);
    Vector<double> chanFreqs(nchan);
    for (uint ch=0; ch<nchan; ++ch) {
      double freq = 100.e6 + ch*10.e6;
      if (ch>2) {
        // Make frequencies unevenly spaced
        freq += 4.e6;
      }
      if (ch>4) {
        freq += 1.e6;
      }
      if (ch>5) {
        freq += 15.e6;
      }
      chanFreqs[ch] = freq;
    }
    if (nchan==2) {
      chanFreqs[0] = 100.e6;
      chanFreqs[1] = 101.e6;
    }
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
    buf.setTime (itsCount*itsTimeInterval + itsFirstTime);
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
  double itsFirstTime;
};



// Class to check result of TestInput run by tests.
class TestOutput: public DPStep
{
public:
  enum tests {WeightsNotChanged=1, DataNotChanged=2, DataChanged=4,
  DataEquals=8, WeightEquals=16};
  TestOutput(int ntime, int nchan, int doTest, bool solsHadFreqAxis=true,
             bool solsHadTimeAxis=true)
    : itsCount(0), itsTimeStep(0), itsNTime(ntime), itsNBl(9), itsNChan(nchan),
      itsNCorr(4), itsTimeInterval(5.), itsDoTest(doTest),
      itsSolsHadFreqAxis(solsHadFreqAxis), itsSolsHadTimeAxis(solsHadTimeAxis)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  { 
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(1,0);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 1.0f, 0.0f);

    vector<double> rightTimes(max(itsNTime, 5));
    rightTimes[0] = 0;
    rightTimes[1] = 2;
    rightTimes[2] = 3;
    for (int t=3; t<itsNTime; ++t) {
      rightTimes[t] = 4;
    }
    if (!itsSolsHadTimeAxis) {
      rightTimes.assign(itsNTime, 0);
    }

    vector<double> rightFreqs(max(itsNChan, 5));
    rightFreqs[0] = 1;
    rightFreqs[1] = 1;
    rightFreqs[2] = 2;
    rightFreqs[3] = 2;
    rightFreqs[4] = 2;
    for (int f=5; f<itsNChan; ++f) {
      rightFreqs[f] = 3;
    }
    if (!itsSolsHadFreqAxis) {
      rightFreqs.assign(itsNChan, 1);
    }

    if (itsDoTest) {
      //cout<<endl;
      for (uint bl=0; bl<info().nbaselines(); ++bl) {
        for (int chan=0; chan<itsNChan; ++chan) {
            uint ant1 = info().getAnt1()[bl];
            uint ant2 = info().getAnt2()[bl];
            // Square root of autocorrelation for first antenna
            complex<float> val = sqrt(buf.getData().data()[bl*itsNCorr*itsNChan + chan*itsNCorr]);

            bool flag = buf.getFlags().data()[bl*itsNCorr*itsNChan + chan*itsNCorr];
            if ((ant1==1 || ant2==1) && rightTimes[itsTimeStep]==2 && rightFreqs[chan]==2) {
              ASSERT(flag);
            } else {
              ASSERT(!flag);
              ASSERT(near(rightTimes[itsTimeStep]*100 + rightFreqs[chan], val));
            }
        }
      }
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
  bool itsSolsHadFreqAxis, itsSolsHadTimeAxis;
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

// Test amplitude correction
void testampl(int ntime, int nchan, bool freqaxis, bool timeaxis)
{
  cout << "testampl: ntime=" << ntime << " nchan=" << nchan << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nchan);
  DPStep::ShPtr step1(in);

  ParameterSet parset1;
  parset1.add ("correction", "myampl");
  parset1.add ("parmdb", "tApplyCalH5_tmp.h5");
  DPStep::ShPtr step2(new ApplyCal(in, parset1, ""));

  DPStep::ShPtr step3(new TestOutput(ntime, nchan,
      TestOutput::WeightsNotChanged, freqaxis, timeaxis));

  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
  cout<<endl;
}


// Write a temporary H5Parm
void createH5Parm(vector<double> times, vector<double> freqs) {
  H5Parm h5parm("tApplyCalH5_tmp.h5", true);

  // Add antenna metadata
  vector<string> antNames;
  vector<vector<double> > antPositions;
  vector<double> oneAntPos(3, 42.);
  for (uint i=0; i<3; ++i) {
    stringstream antNameStr;
    antNameStr<<"ant"<<(i+1);
    antNames.push_back(antNameStr.str());
    antPositions.push_back(oneAntPos);
  }
  h5parm.addAntennas(antNames, antPositions);

  vector<H5Parm::AxisInfo> axes;
  axes.push_back(H5Parm::AxisInfo("ant",3));
  if (!times.empty()) {
    axes.push_back(H5Parm::AxisInfo("time", times.size()));
  }
  if (!freqs.empty()) {
    axes.push_back(H5Parm::AxisInfo("freq", freqs.size()));
  }

  H5Parm::SolTab soltab = h5parm.createSolTab("myampl","amplitude",axes);
  ASSERT(h5parm.nSolTabs() == 1);
  ASSERT(h5parm.hasSolTab("myampl"));
  soltab.setTimes(times);
  soltab.setFreqs(freqs);
  soltab.setAntennas(antNames);

  uint ntimes = max(times.size(), 1);
  uint nfreqs = max(freqs.size(), 1);
  vector<double> values(ntimes*nfreqs*3);
  vector<double> weights(ntimes*nfreqs*3);
  for (uint ant=0; ant<3; ++ant) {
    for (uint t=0; t<ntimes; ++t) {
      for (uint f=0; f<nfreqs; ++f) {
        values[ant*ntimes*nfreqs+t*nfreqs + f] = 1./(100.*(t%100)+(1+f));
        weights[ant*ntimes*nfreqs+t*nfreqs + f] = 1.;
        if (ant==1 && t==2 && f==1) {
          weights[ant*ntimes*nfreqs+t*nfreqs + f] = 0.;
        }
      }
    }
  }
  soltab.setValues(values, weights, "CREATE with DPPP tApplyCalH5");
}

int main()
{
  INIT_LOGGER ("tApplyCalH5");

  vector<double> times;
  times.push_back(4472025742.0);
  times.push_back(4472025745.0);
  times.push_back(4472025747.5);
  times.push_back(4472025748.0);
  times.push_back(4472025762.0);
  vector<double> freqs;
  freqs.push_back(90.e6);
  freqs.push_back(139.e6);
  freqs.push_back(170.e6);

  try {
    createH5Parm(times, freqs);
    testampl(5, 7, true, true);
    createH5Parm(times, freqs);
    testampl(5, 2, true, true);
    createH5Parm(times, vector<double>());
    testampl(8, 9, false, true);
    createH5Parm(vector<double>(), freqs);
    testampl(13, 3, true, false);
    createH5Parm(vector<double>(), vector<double>());
    testampl(9, 2, false, false);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
   return 1;
  }
  return 0;
}
