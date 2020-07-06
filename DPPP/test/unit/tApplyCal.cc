// tApplyCal.cc: Test program for class AORFlagger
// Copyright (C) 2013
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
// $Id: tApplyCal.cc 24221 2013-08-02 12:24:48Z tammo $
//
// @author Tammo Jan Dijkema

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayIO.h>

#include <boost/test/unit_test.hpp>

#include "../../ApplyCal.h"
#include "../../DPInput.h"
#include "../../DPBuffer.h"
#include "../../DPInfo.h"
#include "../../../Common/ParameterSet.h"
#include "../../../Common/StringUtil.h"
#include "../../../Common/StreamUtil.h"

using namespace casacore;
using std::vector;
using DP3::ParameterSet;
using DP3::DPPP::DPInput;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::DPInfo;
using DP3::DPPP::ApplyCal;
using DP3::DPPP::DPStep;

BOOST_AUTO_TEST_SUITE(applycal)

// Simple class to generate input arrays.
// 9 baselines, 3 antennas, 4 correlations
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nchan)
    : itsCount(0), itsNTime(ntime), itsNChan(nchan), itsNBl(9), itsNCorr(4),
      itsTimeInterval(5.)
  {
    info().init (itsNCorr, 0, nchan, ntime, 4472025740.0, itsTimeInterval,
        string(), string());
    // Fill the baseline stations; use 3 stations.
    // So they are called 00 01 02 10 11 12 20 21 22, etc.

    vector<int> ant1(itsNBl);
    vector<int> ant2(itsNBl);
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
    vector<string> antNames {"ant1", "ant2", "ant3", ""};
    // Define their positions (more or less WSRT RT0-3).
    vector<casacore::MPosition> antPos(4);
    vector<double> vals(3);
    vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
    antPos[0] = casacore::MPosition(casacore::Quantum<casacore::Vector<double> >(vals,"m"),
                          casacore::MPosition::ITRF);
    vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
    antPos[1] = casacore::MPosition(casacore::Quantum<casacore::Vector<double> >(vals,"m"),
                          casacore::MPosition::ITRF);
    vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
    antPos[2] = casacore::MPosition(casacore::Quantum<casacore::Vector<double> >(vals,"m"),
                          casacore::MPosition::ITRF);
    vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
    antPos[3] = casacore::MPosition(casacore::Quantum<casacore::Vector<double> >(vals,"m"),
                          casacore::MPosition::ITRF);
    vector<double> antDiam(4, 70.);
    info().set (antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    vector<double> chanWidth(nchan, 1000000.);
    casacore::Vector<double> chanFreqs(nchan);
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
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(1,0);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    weights=1.;

    casacore::Matrix<double> uvw(3, itsNBl);
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
    casacore::Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags (flags);
    casacore::Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
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
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);

    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(1,0);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 1.0f, 0.0f);

    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }

    // The same gain corrections as in tApplyCal_tmp.parmdb
    vector<casacore::Cube<casacore::Complex> > gains(4); // cube for every corr
    for (int corr=0; corr<4; ++corr) {
      gains[corr].resize(casacore::IPosition(3,2,2,3)); //freq,time,ant;
    }

    gains[0](casacore::Slice(0,2),casacore::Slice(0,2),casacore::Slice(0,3))=1;
    gains[1](casacore::Slice(0,2),casacore::Slice(0,2),casacore::Slice(0,3))=0;
    gains[2](casacore::Slice(0,2),casacore::Slice(0,2),casacore::Slice(0,3))=0;
    gains[3](casacore::Slice(0,2),casacore::Slice(0,2),casacore::Slice(0,3))=1;
    // ant2
    gains[0](0,0,1)=2;
    gains[3](0,0,1)=3;
    gains[0](1,1,1)=casacore::Complex(3.,4.);
    // ant3
    gains[2](1,0,2)=.5;

    if (itsDoTest & DataEquals) {
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
      BOOST_CHECK ( casacore::near(buf.getWeights()(0,0,1),4.));
      BOOST_CHECK ( casacore::near(buf.getWeights()(1,0,1),9.));
      BOOST_CHECK ( casacore::near(buf.getWeights()(2,0,1),4.));
      BOOST_CHECK ( casacore::near(buf.getWeights()(3,0,1),9.));
      BOOST_CHECK ( casacore::near(buf.getWeights()(0,31,5),0.8));
    }

    if (itsDoTest & DataEquals) {
      BOOST_CHECK (allNear (buf.getData(), data, 1.e-7));
    }

    if (itsDoTest & DataNotChanged) {
      BOOST_CHECK (allNear (buf.getData(), data, 1.e-7));
    }
    if (itsDoTest & DataChanged) {
      BOOST_CHECK (!(allNear (buf.getData(), data, 1.e-7)));
    }
    if (itsDoTest & WeightsNotChanged) {
      BOOST_CHECK (allNear (buf.getWeights(), weights, 1.e-7));
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
    BOOST_CHECK_EQUAL (itsNChan, int(infoIn.origNChan()));
    BOOST_CHECK_EQUAL (itsNChan, int(infoIn.nchan()));
    BOOST_CHECK_EQUAL (itsNTime, int(infoIn.ntime()));
    BOOST_CHECK_EQUAL (itsTimeInterval, infoIn.timeInterval());
    BOOST_CHECK_EQUAL (itsNBl, int(infoIn.nbaselines()));
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

  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test clock + tec, and test two ApplyCals in sequence
void testclocktec(int ntime, int nchan)
{
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
}

// Test gain
void testgain(int ntime, int nchan)
{
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

BOOST_AUTO_TEST_CASE( test_clock_and_tec ) {
  testclocktec (10,  32);
}

BOOST_AUTO_TEST_CASE( test_gain ) {
  testgain (10, 32);
}

BOOST_AUTO_TEST_SUITE_END()
