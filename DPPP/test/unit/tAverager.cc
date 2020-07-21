// tAverager.cc: Test program for class Averager
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
// $Id$
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayIO.h>

#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>

#include "../../Averager.h"
#include "../../DPBuffer.h"
#include "../../DPInfo.h"
#include "../../../Common/ParameterSet.h"
#include "../../../Common/StringUtil.h"

using std::vector;
using DP3::ParameterSet;
using DP3::DPPP::DPInput;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::DPInfo;
using DP3::DPPP::Averager;
using DP3::DPPP::DPStep;

BOOST_AUTO_TEST_SUITE(averager)

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
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(i+itsCount*10,i-1000+itsCount*6);
    }
    DPBuffer buf;
    buf.setTime (itsCount*5 + 2);   //same interval as in updateAveragInfo
    buf.setData (data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights (weights);
    casacore::Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags (flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    casacore::Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags (fullResFlags);
    casacore::Matrix<double> uvw(3,itsNBl);
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
    info().init (itsNCorr, 0, itsNChan, itsNTime, 100, 5, string(), string());
    // Define the frequencies.
    casacore::Vector<double> chanFreqs(itsNChan);
    vector<double> chanWidth(itsNChan, 100000.);
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
    casacore::Cube<casacore::Complex>  data(itsNCorr,itsNChan,itsNBl);
    casacore::Cube<float> weights(itsNCorr,itsNChan,itsNBl);
    casacore::Cube<bool> fullResFlags(itsNChan,itsNAvgTime,itsNBl);
    fullResFlags = true;   // takes care of missing times at the end
    weights = 0;
    for (int j=itsCount*itsNAvgTime; j<itsCount*itsNAvgTime+navgtime; ++j) {
      for (int i=0; i<int(data.size()); ++i) {
        data.data()[i] += casacore::Complex(i+j*10,i-1000+j*6);
        weights.data()[i] += float(1);
      }
      fullResFlags(casacore::Slicer(casacore::IPosition(3,0,0,0),
                          casacore::IPosition(3,itsNChan,navgtime,itsNBl))) = itsFlag;
    }
    casacore::Cube<casacore::Complex> result(itsNCorr,nchan,itsNBl);
    casacore::Cube<float> resultw(itsNCorr,nchan,itsNBl);
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
    BOOST_CHECK (allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK (allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK (allEQ(buf.getFlags(), itsFlag));
    BOOST_CHECK (casacore::near(buf.getTime(),
                 2+5*(itsCount*itsNAvgTime + (itsNAvgTime-1)/2.)));
    BOOST_CHECK (allNear(buf.getWeights(), resultw, 1e-5));
    if (navgtime == itsNAvgTime) {
      casacore::Matrix<double> uvw(3,itsNBl);
      indgen (uvw, 100*(itsCount*itsNAvgTime + 0.5*(itsNAvgTime-1)));
      BOOST_CHECK (allNear(buf.getUVW(), uvw, 1e-5));
    }
    BOOST_CHECK (allEQ(buf.getFullResFlags(), fullResFlags));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    BOOST_CHECK_EQUAL (itsNChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL (1+(itsNChan-1)/itsNAvgChan, int(info.nchan()));
    BOOST_CHECK_EQUAL (1+(itsNTime-1)/itsNAvgTime, int(info.ntime()));
    BOOST_CHECK_EQUAL (5*itsNAvgTime, info.timeInterval());
    BOOST_CHECK_EQUAL (itsNAvgChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL (itsNAvgTime, int(info.ntimeAvg()));
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
    casacore::Cube<casacore::Complex> data(itsNrCorr,itsNrChan,itsNrBl);
    casacore::Cube<float> weights(itsNrCorr,itsNrChan,itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr,itsNrChan,itsNrBl);
    int i = 0;
    for (int ib=0; ib<itsNrBl; ++ib) {
      for (int ic=0; ic<itsNrChan; ++ic) {
        for (int ip=0; ip<itsNrCorr; ++ip) {
          data.data()[i] = casacore::Complex(i+itsCount*10,i-1000+itsCount*6);
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
    vector<DP3::rownr_t> rownrs(1,itsCount);
    buf.setRowNrs (rownrs);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void getUVW (const casacore::RefRows&, double, DPBuffer& buf)
  {
    buf.getUVW().resize (3, itsNrBl);
    indgen (buf.getUVW());
  }
  virtual bool getFullResFlags (const casacore::RefRows&, DPBuffer& buf)
  {
    buf.getFullResFlags().assign (itsFullResFlags);
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo&)
  {
    // Use timeInterval=5
    info().init (itsNrCorr, 0, itsNrChan, itsNrTime, 100, 5, string(), string());
    // Define the frequencies.
    casacore::Vector<double> chanFreqs(itsNrChan);
    vector<double> chanWidth(itsNrChan, 100000.);
    indgen (chanFreqs, 1050000., 100000.);
    info().set (chanFreqs, chanWidth);
  }
  int itsCount, itsNrTime, itsNrBl, itsNrChan, itsNrCorr;
  casacore::Cube<bool> itsFullResFlags;
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
    casacore::Cube<casacore::Complex> result(itsNrCorr,1,itsNrBl);
    casacore::Cube<float> weights(itsNrCorr,1,itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr,1,itsNrBl);
    casacore::Cube<bool> fullResFlags(itsNrChan,itsNrTime,itsNrBl);
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
              result(ip,0,ib) += weight * casacore::Complex(i+it*10,i-1000+it*6);
              weights(ip,0,ib) += weight;
              flags(ip,0,ib) = false;
              fullResFlags(ic,it,ib) = false;
            }
            i++;
          }
        }
      }
    }
    BOOST_CHECK (allNE(weights, float(0.)));
    for (unsigned int i=0; i<result.size(); ++i) {
      result.data()[i] /= weights.data()[i];
    }
    // Check the averaged result.
    BOOST_CHECK (allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK (allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK (allEQ(buf.getFlags(), flags));
    BOOST_CHECK (casacore::near(buf.getTime(), 2.+5*(itsNrTime-1)/2.));
    BOOST_CHECK (allNear(buf.getWeights(), weights, 1e-5));
    casacore::Matrix<double> uvw(3,itsNrBl);
    indgen (uvw);
    BOOST_CHECK (allNear(buf.getUVW(), uvw, 1e-5));
    BOOST_CHECK (allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    BOOST_CHECK_EQUAL (itsNrChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL (size_t {1}, info.nchan());
    BOOST_CHECK_EQUAL (size_t {1}, info.ntime());
    BOOST_CHECK_EQUAL (5*itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL (itsNrChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL (itsNrTime, int(info.ntimeAvg()));
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
    casacore::Cube<casacore::Complex> result(itsNrCorr,1,itsNrBl);
    casacore::Cube<float> weights(itsNrCorr,1,itsNrBl);
    casacore::Cube<bool> flags(itsNrCorr,1,itsNrBl);
    casacore::Cube<bool> fullResFlags(itsNrChan,itsNrTime,itsNrBl);
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
            i += itsNrCorr;
          } else {
            for (int ip=0; ip<itsNrCorr; ++ip) {
              if ((it+2*ib+3*ic) % 7 != 0) {
                float weight = (1 + (it+ib+ic)%5) / 5.;
                result(ip,0,ib) += weight * casacore::Complex(i+it*10,i-1000+it*6);
                weights(ip,0,ib) += weight;
                flags(ip,0,ib) = false;
                fullResFlags(ic,it,ib) = false;
              }
              i++;
            }
          }
        }
      }
    }
    for (unsigned int i=0; i<result.size(); ++i) {
      if (!flags.data()[i]) {
        result.data()[i] /= weights.data()[i];
      }
    }
    // Check the averaged result.
    BOOST_CHECK (allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK (allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK (allEQ(buf.getFlags(), flags));
    BOOST_CHECK (casacore::near(buf.getTime(), 2.+5*(itsNrTime-1)/2.));
    BOOST_CHECK (allNear(buf.getWeights(), weights, 1e-5));
    casacore::Matrix<double> uvw(3,itsNrBl);
    indgen (uvw);
    BOOST_CHECK (allNear(buf.getUVW(), uvw, 1e-5));
    BOOST_CHECK (allEQ(buf.getFullResFlags(), fullResFlags));
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& info)
  {
    BOOST_CHECK_EQUAL (itsNrChan, int(info.origNChan()));
    BOOST_CHECK_EQUAL (size_t {1}, info.nchan());
    BOOST_CHECK_EQUAL (size_t {1}, info.ntime());
    BOOST_CHECK_EQUAL (5*itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL (itsNrChan, int(info.nchanAvg()));
    BOOST_CHECK_EQUAL (itsNrTime, int(info.ntimeAvg()));
  }

  int itsNrTime, itsNrBl, itsNrChan, itsNrCorr, itsStep;
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

// Test simple averaging without flagged points.
void test1(int ntime, int nbl, int nchan, int ncorr,
           int navgtime, int navgchan, bool flag)
{
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqstep", DP3::toString(navgchan));
  parset.add ("timestep", DP3::toString(navgtime));
  DPStep::ShPtr step2(new Averager(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     navgtime, navgchan, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Like test 1, but specify target resolution
void test1resolution(int ntime, int nbl, int nchan, int ncorr,
                     double timeresolution, double freqresolution,
                     string frequnit, bool flag)
{
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqresolution", DP3::toString(freqresolution)+frequnit);
  parset.add ("timeresolution", DP3::toString(timeresolution));
  DPStep::ShPtr step2(new Averager(in, parset, ""));

  if (!frequnit.empty()) {
    casacore::Quantity q(freqresolution, frequnit);
    freqresolution = q.getValue("Hz", true);
  }

  int navgchan = std::max(1, int(freqresolution / 100000 + 0.5));
  int navgtime = std::max(1, int(timeresolution / 5. + 0.5));

  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     navgtime, navgchan, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Like test1, but the averaging is done in two steps.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag)
{
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

BOOST_AUTO_TEST_CASE( testaverager1 ) {
  test1(10, 3, 32, 4, 2, 4, false);
}

BOOST_AUTO_TEST_CASE( testaverager2 ) {
  test1(10, 3, 30, 1, 3, 3, true);
}

BOOST_AUTO_TEST_CASE( testaverager3 ) {
  test1(10, 3, 30, 1, 3, 3, false);
}
BOOST_AUTO_TEST_CASE( testaverager4 ) {
  test1(11, 3, 30, 2, 3, 3, false);
}

BOOST_AUTO_TEST_CASE( testaverager5 ) {
  test1(10, 3, 32, 4, 1, 32, false);
}

BOOST_AUTO_TEST_CASE( testaverager6 ) {
  test1(10, 3, 32, 1, 1, 1, false);
}

BOOST_AUTO_TEST_CASE( testaverager7 ) {
  test2(10, 3, 32, 2, true);
}

BOOST_AUTO_TEST_CASE( testaverager8 ) {
  test2(10, 3, 32, 2, false);
}

BOOST_AUTO_TEST_CASE( testaverager9 ) {
  test3(1, 1);
}

BOOST_AUTO_TEST_CASE( testaverager10 ) {
  test3(10, 4);
}

BOOST_AUTO_TEST_CASE( testaverager11 ) {
  test4(1, 4, 3);
}

BOOST_AUTO_TEST_CASE( testaverager12 ) {
  test4(20, 4, 5);
}

BOOST_AUTO_TEST_CASE( testresolution1 ) {
  test1resolution(10, 3, 32, 4, 10., 100000, "Hz", false);
}

BOOST_AUTO_TEST_CASE( testresolution2 ) {
  test1resolution(11, 3, 32, 4, 1., 800, "kHz", false);
}

BOOST_AUTO_TEST_CASE( testresolution3 ) {
  test1resolution(11, 3, 32, 4, 15., 0.4, "MHz", false);
}

BOOST_AUTO_TEST_SUITE_END()