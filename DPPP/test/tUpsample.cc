//# tAverager.cc: Test program for class Averager
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
//# $Id: tAverager.cc 35179 2016-08-25 11:25:17Z dijkema $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/Upsample.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>

#include <casa/Quanta/Quantum.h>
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
  TestInput(vector<double> times, vector<bool> flags, double timeInterval)
    : itsTimeStep(0), itsNBl(3), itsNChan(5), itsNCorr(4), itsTimes(times),
      itsFlags(flags), itsTimeInterval(timeInterval)
  {}
private:
  virtual bool process (const DPBuffer&)
  {
    // Stop when all times are done.
    if (itsTimeStep == itsTimes.size()) {
      return false;
    }
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsTimeStep*10,i-1000+itsTimeStep*6);
    }
    DPBuffer buf;
    buf.setTime (itsTimes[itsTimeStep]);
    buf.setData (data);
    Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights (weights);
    Cube<bool> flags(data.shape());
    flags = itsFlags[itsTimeStep];
    buf.setFlags (flags);
    buf.setExposure(itsTimeInterval);

    Matrix<double> uvw(3,itsNBl);
    indgen (uvw, double(itsTimeStep*100));
    buf.setUVW (uvw);
    getNextStep()->process (buf);
    ++itsTimeStep;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}

  virtual void show (std::ostream&) const {}

  virtual void updateInfo (const DPInfo&)
  {
    info().init (itsNCorr, itsNChan, itsTimes.size(), itsTimes[0], itsTimeInterval, string(), string());
    // Define the frequencies.
    Vector<double> chanFreqs(itsNChan);
    Vector<double> chanWidth(itsNChan, 100000.);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
  uint itsTimeStep, itsNBl, itsNChan, itsNCorr;
  vector<double> itsTimes;
  vector<bool> itsFlags;
  double itsTimeInterval;
};

// Class to check result of upsampling TestInput
class TestOutput: public DPStep
{
public:
  TestOutput(vector<double> times, vector<bool> flags, double timeInterval)
    : itsTimes(times), itsFlags(flags), itsTimeStep(0), itsTimeInterval(timeInterval)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    ASSERT(nearAbs(buf.getTime(), itsTimes[itsTimeStep], itsTimeInterval*0.01));
    ASSERT(allTrue(buf.getFlags()) == itsFlags[itsTimeStep]);
    ++itsTimeStep;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    ASSERT(near(info.timeInterval(), itsTimeInterval));
  }

  vector<double> itsTimes;
  vector<bool> itsFlags;
  uint itsTimeStep;
  double itsTimeInterval;
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

void test()
{
  {
    // Create the steps.
    double times_array[] = {5020763030.74, 5020763032.75, 5020763034.76, 5020763035.77, 5020763037.78, 5020763039.8};
    vector<double> times(times_array, times_array+6);
    bool flags_array[] = {false, false, true, false, false, false};
    vector<bool> flags(flags_array, flags_array+6);
    TestInput* in = new TestInput(times, flags, 2.01327);
    DPStep::ShPtr in_step(in);
    ParameterSet parset;
    parset.add ("timestep", "2");

    DPStep::ShPtr upsample(new Upsample(in, parset, ""));

    double newtimes_array[] = {5020763030.23, 5020763031.24, 5020763032.25, 5020763033.25, 5020763034.26, 5020763035.27, 5020763036.27, 5020763037.28, 5020763038.29, 5020763039.29, 5020763040.3};
    vector<double> newtimes(newtimes_array, newtimes_array+11);
    bool newflags_array[] = {false, false, false, false, true, false, false, false, false, false, false};
    vector<bool> newflags(newflags_array, newflags_array+11);

    DPStep::ShPtr out_step(new TestOutput(newtimes, newflags, 0.5 * 2.01327));
    in_step->setNextStep (upsample);
    upsample->setNextStep (out_step);
    execute (in_step);
  }
}


int main()
{
  try {
    test();
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
