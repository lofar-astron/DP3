// tDemixer.cc: Test program for class Demixer
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

#include <boost/test/unit_test.hpp>

#include "../../Demixer.h"
#include "../../DPBuffer.h"
#include "../../DPInfo.h"
#include "../../../Common/ParameterSet.h"
#include "../../../Common/StringUtil.h"

using std::vector;
using DP3::ParameterSet;
using DP3::DPPP::DPInput;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::DPInfo;
using DP3::DPPP::Demixer;
using DP3::DPPP::DPStep;

BOOST_AUTO_TEST_SUITE(demixer)

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
  virtual void updateInfo (DPInfo& info)
    // Use startchan=8 and timeInterval=5
    { info.init (itsNCorr, 8, itsNChan, 0, itsNTime, 5, string(), string()); }

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
    if (!itsFlag) {
      for (int j=itsCount*itsNAvgTime; j<itsCount*itsNAvgTime+navgtime; ++j) {
        for (int i=0; i<int(data.size()); ++i) {
          data.data()[i] += casacore::Complex(i+j*10,i-1000+j*6);
          weights.data()[i] += float(1);
        }
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
  virtual void updateInfo (DPInfo& info)
  {
    BOOST_CHECK_EQUAL (size_t {8}, info.startchan());
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


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set DPInfo.
  DPInfo info;
  DPStep::ShPtr step = step1;
  while (step) {
    step->setInfo (info);
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
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("freqstep", DP3::toString(navgchan));
  parset.add ("timestep", DP3::toString(navgtime));
  parset.add ("sources" , "CasA");
  DPStep::ShPtr step2(new Demixer(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr,
                                     navgtime, navgchan, flag));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

BOOST_AUTO_TEST_CASE( test_demixer1 ) {
  test1(10, 3, 32, 4, 2, 4, false);
}

BOOST_AUTO_TEST_CASE( test_demixer2 ) {
  test1(10, 3, 30, 1, 3, 3, true);
}

BOOST_AUTO_TEST_CASE( test_demixer3 ) {
  test1(10, 3, 30, 1, 3, 3, false);
}

BOOST_AUTO_TEST_CASE( test_demixer4 ) {
  test1(11, 3, 30, 2, 3, 3, false);
}

BOOST_AUTO_TEST_CASE( test_demixer5 ) {
  test1(10, 3, 32, 4, 1, 32, false);
}

BOOST_AUTO_TEST_SUITE_END()
