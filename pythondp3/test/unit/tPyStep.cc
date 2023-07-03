// tPyStep.cc: Test program for the python Step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../PyStep.h"

#include <sstream>

#include <boost/test/unit_test.hpp>

#include <xtensor/xtensor.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../../../steps/test/unit/tStepCommon.h"
#include "../../../steps/test/unit/mock/ThrowStep.h"
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
const std::array<std::size_t, 3> kShape{kNBaselines, kNChannels,
                                        kNCorrelations};
}  // namespace

namespace dp3 {
namespace pythondp3 {

BOOST_AUTO_TEST_SUITE(pystep)

// Simple class to generate input arrays. Weights are 1.
class TestInput final : public steps::MockInput {
 public:
  TestInput() : count_(0) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (count_ == kNTimes) {
      return false;
    }
    buffer->setTime(count_ * 5 + 2);
    buffer->ResizeData(kShape);
    for (int i = 0; i < static_cast<int>(buffer->GetData().size()); ++i) {
      buffer->GetData().data()[i] =
          std::complex<float>(i + count_ * 10, i - count_ * 10);
    }
    buffer->ResizeWeights(kShape);
    buffer->GetWeights().fill(1.0f);
    buffer->ResizeUvw(kNBaselines);
    for (std::size_t i = 0; i < buffer->GetUvw().size(); ++i) {
      buffer->GetUvw().data()[i] = count_ * 100 + i;
    }
    getNextStep()->process(std::move(buffer));
    ++count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }

  void updateInfo(const DPInfo&) override {
    // Use timeInterval=5
    info() = DPInfo(kNCorrelations, kNChannels);
    info().setTimes(100.0, 100.0 + (kNTimes - 1) * 5.0, 5.0);
    // Define the frequencies.
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(kNChannels, 100000.);
    for (int i = 0; i < kNChannels; i++) {
      chan_freqs.push_back(1050000. + i * 100000.);
    }
    info().setChannels(std::move(chan_freqs), std::move(chan_width));
  }
  int count_;
};

// Class to check result of the python step.
class TestOutput final : public dp3::steps::test::ThrowStep {
 public:
  TestOutput() : count_(0) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Fill expected result in similar way as TestInput, but multiplied with
    // factor 2.
    xt::xtensor<std::complex<float>, 3> ref_data(kShape);
    for (int i = 0; i < static_cast<int>(ref_data.size()); ++i) {
      ref_data.data()[i] =
          std::complex<float>(2 * (i + count_ * 10), 2 * (i - count_ * 10));
    }

    // Check whether the "visibility" data in the buffer is indeed
    // multiplied by a factor 2
    BOOST_CHECK(xt::allclose(buffer->GetData(), ref_data));

    // Check whether the weights are divided by 2
    const xt::xtensor<float, 3> ref_weights(kShape, 0.5f);
    BOOST_CHECK(xt::allclose(buffer->GetWeights(), ref_weights));

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
  std::weak_ptr<Step> py_step_weak_ptr;

  // Create the steps.
  {
    auto in_step = std::make_shared<TestInput>();
    // Requires mockpystep to be on the PYTHONPATH!
    ParameterSet parset;
    parset.add("mock.python.module", "mockpystep");
    parset.add("mock.python.class", "MockPyStep");
    parset.add("mock.datafactor", "2");
    parset.add("mock.weightsfactor", "0.5");
    {
      std::shared_ptr<Step> py_step = PyStep::create_instance(parset, "mock.");
      auto out_step = std::make_shared<TestOutput>();

      // Monitor lifetime of the python step
      py_step_weak_ptr = py_step;

      // Check whether print statements in show() method are
      // indeed redirected to output stream
      std::ostringstream output_stream_step;
      py_step->show(output_stream_step);
      BOOST_CHECK_EQUAL(
          output_stream_step.str(),
          "\nMockPyStep\n  data factor:    2.0\n  weights factor: 0.5\n");

      common::Fields provided = py_step->getProvidedFields();
      common::Fields required = py_step->getRequiredFields();
      BOOST_CHECK(provided.Data() && provided.Flags() && provided.Weights());
      BOOST_CHECK(!provided.Uvw());
      BOOST_CHECK(required.Data() && required.Weights());
      BOOST_CHECK((!required.Flags()) && (!required.Uvw()));

      dp3::steps::test::Execute({in_step, py_step, out_step});
    }
    // py_step went out of scope here, but is still reachable following
    // the chain of getNextStep() calls from in_step, and thus should
    // still be alive
    BOOST_TEST(!py_step_weak_ptr.expired());
  }
  // in_step went out of scope here
  // py_step should also no longer be alive now
  BOOST_TEST(py_step_weak_ptr.expired());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace pythondp3
}  // namespace dp3
