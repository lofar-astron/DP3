// tScaleData.cc: Test program for class ScaleData
// Copyright (C) 2020
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
// @author Lars Krombeen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../ScaleData.h"
#include "../DPInput.h"
#include "../DPBuffer.h"
#include "../DPInfo.h"
#include "../../Common/ParameterSet.h"
#include "../../Common/StringUtil.h"
#include "../../Common/StreamUtil.h"

using namespace DP3;
using namespace DP3::DPPP;
using namespace casacore;
using namespace std;

BOOST_AUTO_TEST_SUITE(rotationconstraint)

// Class to check result of TestInput run by test1.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr) {}
public:
  std::unique_ptr<BDABuffer> itsResults;
private:
  virtual bool process (const DPBuffer&) { return true; }
  virtual bool process (std::unique_ptr<BDABuffer> results) { 
    itsResults = std::move(results);
    return true; 
  }
  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn) {}
private:
  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// TODO test with an output step that does not support BDA buffer

BOOST_AUTO_TEST_CASE( test_processing_for_bda_buffer )
{
  int ntime {10};
  int nbl {4};
  int nchan {4};
  int ncorr {4};

  // Preparation
  ParameterSet parset;
  parset.add ("stations", "[rs01.s01, *]");
  parset.add ("coeffs", "[[2,0.5],[3,2,1]]");
  parset.add ("scalesize", "false");
  
  DPInfo info = DPInfo();
  info.init (ncorr, 0, nchan, ntime, 0., 5., string(), string());
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
  antPos[0] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
  antPos[1] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
  antPos[2] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
  antPos[3] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  Vector<double> antDiam(4, 70.);
  info.set(antNames, antDiam, antPos, ant1, ant2);
  // Define the frequencies.
  Vector<double> chanWidth(nchan, 1000000.);
  Vector<double> chanFreqs(nchan);
  indgen (chanFreqs, 10500000., 1000000.);
  info.set (chanFreqs, chanWidth);

  std::shared_ptr<ScaleData> stepScaleData(new ScaleData(nullptr, parset, ""));
  std::shared_ptr<TestOutput> stepTestOutput(new TestOutput(ntime, nbl, nchan, ncorr));
  stepScaleData->setNextStep (stepTestOutput);
  stepScaleData->setInfo(info);

  // Initialize buffer
  const int datasize {nbl * nchan * ncorr};
  std::unique_ptr<BDABuffer> bdaBuffer { new BDABuffer(datasize) };
  for (int ind = 0; ind < nbl; ++ind)
  {
    const std::complex<float> data = ind + 2;
    bdaBuffer->addRow(ntime, 5., ind, nchan, ncorr, &data, nullptr, nullptr, nullptr, nullptr);
  }

  // // Execution
  stepScaleData->process(std::move(bdaBuffer));

  // Assertion
  const auto results = stepTestOutput->itsResults->getData();
  for (int ind = 0; ind < datasize; ++ind)
  {
    const std::complex<float> data = ind + 2;
    BOOST_CHECK_EQUAL(data * 2, results[ind]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
