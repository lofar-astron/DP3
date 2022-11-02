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
#include "mock/ThrowStep.h"
#include "../../../pythondp3/PyStep.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Step;

namespace {
const int kNTimes = 10;
const int kNBaselines = 3;
const int kNChannels = 32;
const int kNCorrelations = 4;
}  // namespace

namespace dp3 {
namespace pythondp3 {

BOOST_AUTO_TEST_SUITE(pystep)

// Simple class to generate input arrays. Weights are 1.
class TestInput final : public steps::MockInput {
 public:
  TestInput() : count_(0) {}

 private:
  bool process(const DPBuffer&) override {
    // Stop when all times are done.
    if (count_ == kNTimes) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(kNCorrelations, kNChannels,
                                           kNBaselines);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(i + count_ * 10, i - count_ * 10);
    }
    DPBuffer buf;
    buf.setTime(count_ * 5 + 2);
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Matrix<double> uvw(3, kNBaselines);
    indgen(uvw, double(count_ * 100));
    buf.setUVW(uvw);
    getNextStep()->process(buf);
    ++count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Use timeInterval=5
    info().init(kNCorrelations, 0, kNChannels, kNTimes, 100, 5, std::string());
    // Define the frequencies.
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(kNChannels, 100000.);
    for (int i = 0; i < kNChannels; i++) {
      chan_freqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chan_freqs), std::move(chan_width));
  }
  int count_;
};

// Class to check result of the python step.
class TestOutput final : public dp3::steps::test::ThrowStep {
 public:
  TestOutput() : count_(0) {}

  bool process(const DPBuffer& buf) override {
    // Fill expected result in similar way as TestInput, but multiplied with
    // factor 2.
    casacore::Cube<casacore::Complex> ref_data(kNCorrelations, kNChannels,
                                               kNBaselines);
    casacore::Cube<float> weights(kNCorrelations, kNChannels, kNBaselines);

    for (int i = 0; i < int(ref_data.size()); ++i) {
      ref_data.data()[i] =
          casacore::Complex(2 * (i + count_ * 10), 2 * (i - count_ * 10));
    }

    casacore::Cube<float> ref_weights(kNCorrelations, kNChannels, kNBaselines);
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
  void updateInfo(const DPInfo&) override {}

 private:
  int count_;
};

BOOST_AUTO_TEST_CASE(simple_pystep) {
  // Weak pointer that will be used to monitor the lifetime of the last step
  // Because of the complicated ownership across the Python-C++ boundary
  // we need to check whether the steps are destroyed at the right time
  std::weak_ptr<Step> out_weak_ptr;

  // Create the steps.
  {
    auto in_step = std::make_shared<TestInput>();
    // Requires mockpystep to be on the PYTHONPATH!
    ParameterSet parset;
    parset.add("python.module", "mockpystep");
    parset.add("python.class", "MockPyStep");
    parset.add("datafactor", "2");
    parset.add("weightsfactor", "0.5");
    {
      std::shared_ptr<Step> py_step =
          PyStep::create_instance(in_step.get(), parset, "");
      auto out_step = std::make_shared<TestOutput>();

      // Monitor lifetime of output step
      out_weak_ptr = out_step;

      // Check whether print statements in show() method are
      // indeed redirected to output stream
      std::ostringstream output_stream_step;
      py_step->show(output_stream_step);
      BOOST_TEST(
          output_stream_step.str() ==
          "\nMockPyStep\n  data factor:    2.0\n  weights factor: 0.5\n");

      dp3::steps::test::Execute({in_step, py_step, out_step});
    }
    // out_step went out of scope here, but is still reachable following
    // the chain of getNextStep() calls from in_step, and thus should
    // still be alive
    BOOST_TEST(!out_weak_ptr.expired());
  }
  // in_step went out of scope here
  // out_step should also no longer be alive now
  BOOST_TEST(out_weak_ptr.expired());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace pythondp3
}  // namespace dp3
