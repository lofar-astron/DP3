// tPyDPStep.cc: Test program for the python DPStep
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>

#include "../../../PythonDPPP/PyDPStep.h"
#include "../../DPBuffer.h"
#include "../../DPInfo.h"
#include "../../DPInput.h"
#include "../../../Common/ParameterSet.h"

using DP3::ParameterSet;
using DP3::DPPP::DPBuffer;
using DP3::DPPP::DPInfo;
using DP3::DPPP::DPInput;
using DP3::DPPP::DPStep;
using DP3::DPPP::PyDPStep;
using std::vector;

BOOST_AUTO_TEST_SUITE(pydpstep)

// Simple class to generate input arrays.
// Weights are 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public DPInput {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
      : count_(0),
        ntimes_(ntime),
        nblines_(nbl),
        nchan_(nchan),
        ncorr_(ncorr) {}

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (count_ == ntimes_) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(ncorr_, nchan_, nblines_);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(i + count_ * 10, i - count_ * 10);
    }
    DPBuffer buf;
    buf.setTime(count_ * 5 + 2);
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Matrix<double> uvw(3, nblines_);
    indgen(uvw, double(count_ * 100));
    buf.setUVW(uvw);
    getNextStep()->process(buf);
    ++count_;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {
    // Use timeInterval=5
    info().init(ncorr_, 0, nchan_, ntimes_, 100, 5, string(), string());
    // Define the frequencies.
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(nchan_, 100000.);
    for (int i = 0; i < nchan_; i++) {
      chan_freqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chan_freqs), std::move(chan_width));
  }
  int count_, ntimes_, nblines_, nchan_, ncorr_;
};

// Class to check result of the python step.
class TestOutput : public DPStep {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
      : count_(0),
        ntimes_(ntime),
        nblines_(nbl),
        nchan_(nchan),
        ncorr_(ncorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Fill expected result in similar way as TestInput, but multiplied with
    // factor 2.
    casacore::Cube<casacore::Complex> ref_data(ncorr_, nchan_, nblines_);
    casacore::Cube<float> weights(ncorr_, nchan_, nblines_);

    for (int i = 0; i < int(ref_data.size()); ++i) {
      ref_data.data()[i] =
          casacore::Complex(2 * (i + count_ * 10), 2 * (i - count_ * 10));
    }

    casacore::Cube<float> ref_weights(ncorr_, nchan_, nblines_);
    ref_weights = 0.5;

    // Check whether the "visibility" data in the buffer is indeed
    // multiplied by a factor 2
    BOOST_CHECK(allNear(real(buf.getData()), real(ref_data), 1e-5));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(ref_data), 1e-5));

    // Check whether the weights are divided by 2
    BOOST_CHECK(allNear(buf.getWeights(), ref_weights, 1e-5));

    ++count_;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& info) {}

  int count_;
  int ntimes_, nblines_, nchan_, ncorr_, itsNAvgTime, itsNAvgChan;
};

// Execute steps.
void Execute(const DPStep::ShPtr& step1) {
  // Set DPInfo.
  step1->setInfo(DPInfo());
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf))
    ;
  step1->finish();
}

// Test simple python step
void test(int ntime, int nbl, int nchan, int ncorr) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  // Requires mockpystep to be on the PYTHONPATH!
  ParameterSet parset;
  parset.add("python.module", "mockpystep");
  parset.add("python.class", "MockPyStep");
  parset.add("datafactor", "2");
  parset.add("weightsfactor", "0.5");
  // Step 2 is the python step
  DPStep::ShPtr step2 = PyDPStep::create_instance(in, parset, "");
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr));
  step1->setNextStep(step2);
  step2->setNextStep(step3);

  // Check whether print statements in show() method are
  // indeed redirected to output stream
  std::ostringstream output_stream_step;
  step2->show(output_stream_step);
  BOOST_TEST(output_stream_step.str() ==
             "\nMockPyStep\n  data factor:    2.0\n  weights factor: 0.5\n");

  Execute(step1);
}

BOOST_AUTO_TEST_CASE(testpydpstep) { test(10, 3, 32, 4); }

BOOST_AUTO_TEST_SUITE_END()