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
    : count_(0), ntime_(ntime), nbl_(nbl), nchan_(nchan),
      ncorr_(ncorr) {}
public:
  std::unique_ptr<BDABuffer> results_;
private:
  virtual bool process (const DPBuffer&) { return true; }
  virtual bool process (std::unique_ptr<BDABuffer> results) {
    results_ = std::move(results);
    return true;
  }
  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn) {}
private:
  int count_;
  int ntime_, nbl_, nchan_, ncorr_;
};

// Generate DP Info for 2 antennas
DPInfo GenerateDPInfo(int ntime, int nbl, int nchan, int ncorr) {
  DPInfo info = DPInfo();
  info.init (ncorr, 0, nchan, ntime, 0., 5., string(), string());
  // Fill the baseline stations; use 2 stations.
  // So they are called 00 01 02 03 10 11 12 13 20, etc.
  Vector<Int> ant1(nbl);
  Vector<Int> ant2(nbl);
  int st1 = 0;
  int st2 = 0;
  for (int i=0; i<nbl; ++i) {
    ant1[i] = st1;
    ant2[i] = st2;
    if (++st2 == 2) {
      st2 = 0;
      if (++st1 == 2) {
        st1 = 0;
      }
    }
  }
  Vector<String> ant_names(2);
  ant_names[0] = "rs01.s01";
  ant_names[1] = "rs02.s01";
  // Define their positions (more or less WSRT RT0-3).
  vector<MPosition> ant_pos(2);
  Vector<double> vals(3);
  vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
  ant_pos[0] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
  ant_pos[1] = MPosition(Quantum<Vector<double> >(vals,"m"), MPosition::ITRF);
  Vector<double> antDiam(2, 70.);
  info.set(ant_names, antDiam, ant_pos, ant1, ant2);
  // Define the frequencies.
  Vector<double> chan_width(nchan, 3.);
  Vector<double> chan_freqs(nchan);
  indgen (chan_freqs, 1., 3.);
  info.set (chan_freqs, chan_width);

  return info;
}

BOOST_AUTO_TEST_CASE( test_processing_for_bda_buffer )
{
  int ntime {10};
  int nbl {2};
  int nchan {1};
  int ncorr {2};
  int nantennas {2};

  // Preparation
  ParameterSet parset;
  parset.add ("stations", "[rs01.s01, *]");
  parset.add ("coeffs", "[[2,1],[3,2,1]]");
  parset.add ("scalesize", "false");

  DPInfo info = GenerateDPInfo(ntime, nbl, nchan, ncorr);

  std::shared_ptr<ScaleData> step_scale_data(new ScaleData(nullptr, parset, ""));
  std::shared_ptr<TestOutput> step_test_output(new TestOutput(ntime, nbl, nchan, ncorr));
  step_scale_data->setNextStep (step_test_output);
  step_scale_data->setInfo(info);

  // Initialize buffer
  const int datasize {nbl * nchan * ncorr};
  std::unique_ptr<BDABuffer> bda_buffer { new BDABuffer(datasize) };
  for (int ind = 0; ind < nbl * nantennas; ++ind)
  {
    const std::complex<float> data = ind + 1;
    bda_buffer->AddRow(ntime, 5., 0, nchan, ncorr, ind % nantennas, &data, nullptr, nullptr, nullptr, nullptr);
  }

  // // Execution
  step_scale_data->process(std::move(bda_buffer));

  // Assertion
  const auto results = step_test_output->results_->GetData();
  // size shoule be equal to datasize
  BOOST_CHECK_EQUAL(size_t{4}, step_test_output->results_->GetNumberOfElements());
  BOOST_CHECK(near(4., results[0].real()));
  BOOST_CHECK(near(9.798, results[2].real()));
  // Results 1 and 3 are close to zero, but slightly different every test.
}

BOOST_AUTO_TEST_SUITE_END()
