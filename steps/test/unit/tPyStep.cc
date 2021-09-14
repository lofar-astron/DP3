// tPyStep.cc: Test program for the python Step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../../pythondp3/PyStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../InputStep.h"
#include "../../../common/ParameterSet.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::Step;
using std::vector;

namespace dp3 {
namespace pythondp3 {

BOOST_AUTO_TEST_SUITE(pystep)

// Simple class to generate input arrays.
// Weights are 1.
// It can be used with different nr of times, channels, etc.
class TestInput final : public InputStep {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
      : count_(0),
        ntimes_(ntime),
        nblines_(nbl),
        nchan_(nchan),
        ncorr_(ncorr) {}

 private:
  bool process(const DPBuffer&) override {
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

  void finish() override { getNextStep()->finish(); }
  void show(std::ostream&) const override {}
  void updateInfo(const DPInfo&) override {
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
class TestOutput final : public Step {
 public:
  TestOutput(int nbl, int nchan, int ncorr)
      : count_(0), nblines_(nbl), nchan_(nchan), ncorr_(ncorr) {}

 private:
  bool process(const DPBuffer& buf) override {
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

  void finish() override {}
  void show(std::ostream&) const override {}
  void updateInfo(const DPInfo&) override {}

  int count_;
  int nblines_, nchan_, ncorr_;
};

// Test simple python step
void test(int ntime, int nbl, int nchan, int ncorr) {
  // Weak pointer that will be used to monitor the lifetime of the last step
  // Because of the complicated ownership across the Python-C++ boundary
  // we need to check whether the steps are destroyed at the right time
  std::weak_ptr<Step> step3_weak_ptr;

  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);

  {
    Step::ShPtr step1(in);
    // Requires mockpystep to be on the PYTHONPATH!
    ParameterSet parset;
    parset.add("python.module", "mockpystep");
    parset.add("python.class", "MockPyStep");
    parset.add("datafactor", "2");
    parset.add("weightsfactor", "0.5");
    {
      // Step 2 is the python step
      Step::ShPtr step2 = PyStep::create_instance(in, parset, "");
      Step::ShPtr step3(new TestOutput(nbl, nchan, ncorr));

      // Monitor lifetime of output step
      step3_weak_ptr = std::weak_ptr<Step>(step3);

      // Check whether print statements in show() method are
      // indeed redirected to output stream
      std::ostringstream output_stream_step;
      step2->show(output_stream_step);
      BOOST_TEST(
          output_stream_step.str() ==
          "\nMockPyStep\n  data factor:    2.0\n  weights factor: 0.5\n");

      dp3::steps::test::Execute({step1, step2, step3});
    }
    // step3 went out of scope here, but is still reachable following
    // the chain of getNextStep() calls from step1, and thus should
    // still be alive
    BOOST_TEST(!step3_weak_ptr.expired());
  }
  // step1 went out of scope here
  // step3 should also no longer be alive now
  BOOST_TEST(step3_weak_ptr.expired());
}

BOOST_AUTO_TEST_CASE(testpystep) { test(10, 3, 32, 4); }

BOOST_AUTO_TEST_SUITE_END()

}  // namespace pythondp3
}  // namespace dp3
