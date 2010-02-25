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
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/Averager.h>
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
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsFlag(flag)
  {}
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
    indgen (uvw, double(itsCount*100));
    buf.setUVW (uvw);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
    // Use startchan=8 and timeInterval=5
    { avgInfo.init (itsNCorr, 8, itsNChan, itsNBl, itsNTime, 5); }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of averaging TestInput.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr,
             int navgtime, int navgchan, bool flag)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr), itsNAvgTime(navgtime), itsNAvgChan(navgchan),
      itsFlag(flag)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    int nchan = 1+(itsNChan-1)/itsNAvgChan;
    int navgtime = std::min(itsNAvgTime, itsNTime-itsCount*itsNAvgTime);
    // Fill expected result in similar way as TestInput.
    Cube<Complex>  data(itsNCorr,itsNChan,itsNBl);
    Cube<float> weights(itsNCorr,itsNChan,itsNBl);
    Cube<bool> fullResFlags(itsNChan,itsNAvgTime,itsNBl);
    fullResFlags = true;   // takes care of missing times at the end
    weights = 0;
    if (!itsFlag) {
      for (int j=itsCount*itsNAvgTime; j<itsCount*itsNAvgTime+navgtime; ++j) {
        for (int i=0; i<int(data.size()); ++i) {
          data.data()[i] += Complex(i+j*10,i-1000+j*6);
          if (!itsFlag) {
            weights.data()[i] += float(1);
          }
        }
      }
      fullResFlags(Slicer(IPosition(3,0,0,0),
                          IPosition(3,itsNChan,navgtime,itsNBl))) = itsFlag;
    }
    Cube<Complex> result(itsNCorr,nchan,itsNBl);
    Cube<float> resultw(itsNCorr,nchan,itsNBl);
    resultw = 0;
    // Average to get the true expected result.
    for (int k=0; k<itsNBl; ++k) {
      for (int i=0; i<itsNCorr; ++i) {
        for (int j=0; j<nchan; ++j) {
          int jc;
          for (jc=j*itsNAvgChan;
               jc<std::min((j+1)*itsNAvgChan, itsNChan); ++jc) {
            result(i,j,k) += data(i,jc,k);
            resultw(i,j,k) += weights(i,jc,k);
          }
          result(i,j,k) /= float(navgtime*(jc-j*itsNAvgChan));
        }
      }
    }
    // Check the averaged result.
    ASSERT (allNear(real(buf.getData()), real(result), 1e-5));
    ///cout << imag(buf.getData()) << endl<<imag(result);
    ASSERT (allNear(imag(buf.getData()), imag(result), 1e-5));
    ASSERT (allEQ(buf.getFlags(), itsFlag));
    ASSERT (near(buf.getTime(),
                 2+5*(itsCount*itsNAvgTime + (itsNAvgTime-1)/2.)));
    ASSERT (allNear(buf.getWeights(), resultw, 1e-5));
    if (navgtime == itsNAvgTime) {
      Matrix<double> uvw(3,itsNBl);
      indgen (uvw, 100*(itsCount*itsNAvgTime + 0.5*(itsNAvgTime-1)));
      ASSERT (allNear(buf.getUVW(), uvw, 1e-5));
    }
    ///cout <<buf.getFullResFlags()<< fullResFlags;
    ASSERT (allEQ(buf.getFullResFlags(), fullResFlags));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==8);
    ASSERT (int(avgInfo.origNChan())==itsNChan);
    ASSERT (int(avgInfo.nchan())==1+(itsNChan-1)/itsNAvgChan);
    ASSERT (int(avgInfo.ntime())==1+(itsNTime-1)/itsNAvgTime);
    ASSERT (avgInfo.timeInterval()==5*itsNAvgTime);
    ASSERT (int(avgInfo.nchanAvg())==itsNAvgChan);
    ASSERT (int(avgInfo.ntimeAvg())==itsNAvgTime);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
  bool itsFlag;
};


// More elaborate class which can set different flags and weights.
class TestInput3: public DPInput
{
public:
  TestInput3(int nrtime, int nrbl, int nrchan, int nrcorr)
    : itsCount(0),
      itsNrTime(nrtime), itsNrBl(nrbl), itsNrChan(nrchan), itsNrCorr(nrcorr)
  {
    itsFullResFlags.resize (itsNrChan,1,nrbl);
  }
private:
  virtual bool process (const DPBuffer&)
  {
    // Stop when all times are done.
    if (itsCount == itsNrTime) {
      return false;
    }
    Cube<Complex> data(itsNrCorr,itsNrChan,itsNrBl);
    Cube<float> weights(itsNrCorr,itsNrChan,itsNrBl);
    Cube<bool> flags(itsNrCorr,itsNrChan,itsNrBl);
    int i = 0;
    for (int ib=0; ib<itsNrBl; ++ib) {
      for (int ic=0; ic<itsNrChan; ++ic) {
        for (int ip=0; ip<itsNrCorr; ++ip) {
          data.data()[i] = Complex(i+itsCount*10,i-1000+itsCount*6);
          weights.data()[i] = (1 + (itsCount+ib+ic)%5) / 5.;
          flags.data()[i] = ((itsCount+2*ib+3*ic) % 7 == 0);
          i++;
        }
        itsFullResFlags(ic,0,ib) = ((itsCount+2*ib+3*ic) % 7 == 0);
      }
    }
    DPBuffer buf;
    buf.setTime (itsCount*5 + 2);   //same interval as in updateAveragInfo
    buf.setData (data);
    buf.setWeights (weights);
    buf.setFlags (flags);
    Vector<uint> rownrs(1,itsCount);
    buf.setRowNrs (rownrs);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual casa::Matrix<double> getUVW (const casa::RefRows&)
  {
    Matrix<double> uvw(3,itsNrBl);
    indgen (uvw);
    return uvw;
  }
  virtual casa::Cube<bool> getFullResFlags (const casa::RefRows&)
  {
    return itsFullResFlags;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
    // Use startchan=0 and timeInterval=5
    { avgInfo.init (itsNrCorr, 0, itsNrChan, itsNrBl, itsNrTime, 5); }

  int itsCount, itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
  Cube<bool> itsFullResFlags;
};

// Class to check result of averaging TestInput3.
// All input must be averaged (in one or more steps) to a single value
// per corr/baseline.
class TestOutput3: public DPStep
{
public:
  TestOutput3(int nrtime, int nrbl, int nrchan, int nrcorr)
    : itsNrTime(nrtime), itsNrBl(nrbl), itsNrChan(nrchan), itsNrCorr(nrcorr)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    Cube<Complex> result(itsNrCorr,1,itsNrBl);
    Cube<float> weights(itsNrCorr,1,itsNrBl);
    Cube<bool> flags(itsNrCorr,1,itsNrBl);
    Cube<bool> fullResFlags(itsNrChan,itsNrTime,itsNrBl);
    weights = float(0);
    flags = true;
    fullResFlags = true;
    // Create data in the same way as in TestInput3.
    for (int it=0; it<itsNrTime; ++it) {
      int i = 0;
      for (int ib=0; ib<itsNrBl; ++ib) {
        for (int ic=0; ic<itsNrChan; ++ic) {
          for (int ip=0; ip<itsNrCorr; ++ip) {
            if ((it+2*ib+3*ic) % 7 != 0) {
              float weight = (1 + (it+ib+ic)%5) / 5.;
              result(ip,0,ib) += weight * Complex(i+it*10,i-1000+it*6);
              weights(ip,0,ib) += weight;
              ///  cout << result(ip,0,ib)  << weight << endl;
              flags(ip,0,ib) = false;
              fullResFlags(ic,it,ib) = false;
            }
            i++;
          }
        }
      }
    }
    for (uint i=0; i<result.size(); ++i) {
      result.data()[i] /= weights.data()[i];
    }
    // Check the averaged result.
    ///cout << real(buf.getData()) << endl<<real(result);
    ASSERT (allNear(real(buf.getData()), real(result), 1e-5));
    ASSERT (allNear(imag(buf.getData()), imag(result), 1e-5));
    ASSERT (allEQ(buf.getFlags(), flags));
    ASSERT (near(buf.getTime(), 2.+5*(itsNrTime-1)/2.));
    ASSERT (allNear(buf.getWeights(), weights, 1e-5));
    Matrix<double> uvw(3,itsNrBl);
    indgen (uvw);
    ASSERT (allNear(buf.getUVW(), uvw, 1e-5));
    ///cout <<buf.getFullResFlags()<< fullResFlags;
    ASSERT (allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==0);
    ASSERT (int(avgInfo.origNChan())==itsNrChan);
    ASSERT (avgInfo.nchan()==1);
    ASSERT (avgInfo.ntime()==1);
    ASSERT (avgInfo.timeInterval()==5*itsNrTime);
    ASSERT (int(avgInfo.nchanAvg())==itsNrChan);
    ASSERT (int(avgInfo.ntimeAvg())==itsNrTime);
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
};

// Simple class to flag every step-th XX point.
class TestFlagger: public DPStep
{
public:
  TestFlagger(int step)
    : itsCount(0), itsStep(step)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    DPBuffer buf2(buf);
    int ncorr = buf2.getFlags().shape()[0];
    int np = buf2.getFlags().size() / ncorr;
    bool* flagPtr = buf2.getFlags().data();
    for (int i=0; i<np; ++i) {
      if ((i+itsCount)%itsStep == 0) {
        ///cout << "flagged " <<itsCount <<' '<<  i << endl;
        for (int j=0; j<ncorr; ++j) {
          flagPtr[i*ncorr + j] = true;
        }
      }
    }
    getNextStep()->process (buf2);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}

  int itsCount, itsStep;
};

// Class to check result of averaging and flagging TestInput3.
// First the data are averaged from 8,4 to 4,2, then every step-th point
// is flagged, and finally it is averaged to 1,1.
class TestOutput4: public DPStep
{
public:
  TestOutput4(int nrtime, int nrbl, int nrchan, int nrcorr, int step)
    : itsNrTime(nrtime), itsNrBl(nrbl), itsNrChan(nrchan), itsNrCorr(nrcorr),
      itsStep(step)
  {}
private:
  virtual bool process (const DPBuffer& buf)
  {
    Cube<Complex> result(itsNrCorr,1,itsNrBl);
    Cube<float> weights(itsNrCorr,1,itsNrBl);
    Cube<bool> flags(itsNrCorr,1,itsNrBl);
    Cube<bool> fullResFlags(itsNrChan,itsNrTime,itsNrBl);
    weights = float(0);
    flags = true;
    fullResFlags = true;
    // Create data in the same way as in TestInput3.
    for (int it=0; it<itsNrTime; ++it) {
      int i = 0;
      for (int ib=0; ib<itsNrBl; ++ib) {
        for (int ic=0; ic<itsNrChan; ++ic) {
          // TestFlagger flags every step-th point of 2x2 averaged data.
          int tf = it/2;    // same as itsCount in testFlagger
          if (((ib*itsNrChan + ic)/2 + tf) % itsStep == 0) {
            ///cout << "out4 flagged "<< tf<<' '<< i/itsNrCorr<<' ' <<ib<<' '<<ic/2 << endl;
            i += itsNrCorr;
          } else {
            for (int ip=0; ip<itsNrCorr; ++ip) {
              if ((it+2*ib+3*ic) % 7 != 0) {
                float weight = (1 + (it+ib+ic)%5) / 5.;
                result(ip,0,ib) += weight * Complex(i+it*10,i-1000+it*6);
                weights(ip,0,ib) += weight;
                ///  cout << result(ip,0,ib)  << weight << endl;
                flags(ip,0,ib) = false;
                fullResFlags(ic,it,ib) = false;
              }
              i++;
            }
          }
        }
      }
    }
    for (uint i=0; i<result.size(); ++i) {
      if (!flags.data()[i]) {
        result.data()[i] /= weights.data()[i];
      }
    }
    // Check the averaged result.
    ///cout << real(buf.getData()) << endl<<real(result);
    ASSERT (allNear(real(buf.getData()), real(result), 1e-5));
    ASSERT (allNear(imag(buf.getData()), imag(result), 1e-5));
    ASSERT (allEQ(buf.getFlags(), flags));
    ASSERT (near(buf.getTime(), 2.+5*(itsNrTime-1)/2.));
    ASSERT (allNear(buf.getWeights(), weights, 1e-5));
    Matrix<double> uvw(3,itsNrBl);
    indgen (uvw);
    ASSERT (allNear(buf.getUVW(), uvw, 1e-5));
    ///cout <<buf.getFullResFlags()<< fullResFlags;
    ASSERT (allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateAverageInfo (AverageInfo& avgInfo)
  {
    ASSERT (avgInfo.startChan()==0);
    ASSERT (int(avgInfo.origNChan())==itsNrChan);
    ASSERT (avgInfo.nchan()==1);
    ASSERT (avgInfo.ntime()==1);
    ASSERT (avgInfo.timeInterval()==5*itsNrTime);
    ASSERT (int(avgInfo.nchanAvg())==itsNrChan);
    ASSERT (int(avgInfo.ntimeAvg())==itsNrTime);
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr, itsStep;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set AverageInfo.
  AverageInfo avgInfo;
  DPStep::ShPtr step = step1;
  while (step) {
    step->updateAverageInfo (avgInfo);
    step = step->getNextStep();
  }
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test simple averaging without flagged points.
void test1(int ntime, int nbl, int nchan, int ncorr,
           int navgtime, int navgchan, bool flag)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " navgtime=" << navgtime
       << " navgchan=" << navgchan << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqstep", toString(navgchan));
  parset.add ("timestep", toString(navgtime));
  DPStep::ShPtr step2(new Averager(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     navgtime, navgchan, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Like test1, but the averaging is done in two steps.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << " navgtime=2"
       << " navgchan=4" << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset1, parset2;
  parset1.add ("freqstep", "4");
  parset2.add ("timestep", "2");
  DPStep::ShPtr step2a(new Averager(in, parset1, ""));
  DPStep::ShPtr step2b(new Averager(in, parset2, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr, 2, 4, flag));
  step1->setNextStep (step2a);
  step2a->setNextStep (step2b);
  step2b->setNextStep (step3);
  execute (step1);
}

// Do tests with weighting and some flagged points.
void test3(int nrbl, int nrcorr)
{
  {
    cout << "test3: ntime=2 nrbl=" << nrbl << " nchan=2 ncorr=" << nrcorr
         << endl;
    cout << "  navgtime=2 navgchan=2" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(2, nrbl, 2, nrcorr);
    DPStep::ShPtr step1(in);
    ParameterSet parset1;
    parset1.add ("freqstep", "2");
    parset1.add ("timestep", "2");
    DPStep::ShPtr step2a(new Averager(in, parset1, ""));
    DPStep::ShPtr step3(new TestOutput3(2, nrbl, 2, nrcorr));
    step1->setNextStep (step2a);
    step2a->setNextStep (step3);
    execute (step1);
  }
  {
    cout << "test3: ntime=4 nrbl=" << nrbl << " nchan=8 ncorr=" << nrcorr
         << endl;
    cout << "  [navgtime=2 navgchan=4], [navgtime=2 navgchan=2]" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(4, nrbl, 8, nrcorr);
    DPStep::ShPtr step1(in);
    ParameterSet parset1, parset2;
    parset1.add ("freqstep", "4");
    parset1.add ("timestep", "2");
    parset2.add ("freqstep", "2");
    parset2.add ("timestep", "2");
    DPStep::ShPtr step2a(new Averager(in, parset1, ""));
    DPStep::ShPtr step2b(new Averager(in, parset2, ""));
    DPStep::ShPtr step3(new TestOutput3(4, nrbl, 8, nrcorr));
    step1->setNextStep (step2a);
    step2a->setNextStep (step2b);
    step2b->setNextStep (step3);
    execute (step1);
  }
}

// Do tests with averaging and flagging steps to see if the flags are
// promoted to the FULLRES flags.
void test4(int nrbl, int nrcorr, int flagstep)
{
  {
    cout << "test4: ntime=4 nrbl=" << nrbl << " nchan=8 ncorr=" << nrcorr
         << endl;
    cout << "  [navgtime=2 navgchan=2], [flagstep=" << flagstep
         << "] [navgtime=2 navgchan=4]" << endl;
    // Create the steps.
    TestInput3* in = new TestInput3(4, nrbl, 8, nrcorr);
    DPStep::ShPtr step1(in);
    ParameterSet parset1, parset2;
    parset1.add ("freqstep", "2");
    parset1.add ("timestep", "2");
    parset2.add ("freqstep", "4");
    parset2.add ("timestep", "2");
    DPStep::ShPtr step2a(new Averager(in, parset1, ""));
    DPStep::ShPtr step2b(new TestFlagger(flagstep));
    DPStep::ShPtr step2c(new Averager(in, parset2, ""));
    DPStep::ShPtr step3(new TestOutput4(4, nrbl, 8, nrcorr, flagstep));
    step1->setNextStep (step2a);
    step2a->setNextStep (step2b);
    step2b->setNextStep (step2c);
    step2c->setNextStep (step3);
    execute (step1);
  }
}


int main()
{
  try {
    test1(10, 3, 32, 4, 2, 4, false);
    test1(10, 3, 30, 1, 3, 3, true);
    test1(10, 3, 30, 1, 3, 3, false);
    test1(11, 3, 30, 2, 3, 3, false);
    test1(10, 3, 32, 4, 1, 32, false);
    test2(10, 3, 32, 2, true);
    test2(10, 3, 32, 2, false);
    test3(1, 1);
    test3(10, 4);
    test4(1, 4, 3);
    test4(20, 4, 5);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
