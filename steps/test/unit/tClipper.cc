// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/Clipper.h"

#include <iostream>
#include <filesystem>
#include <fstream>

#include <boost/test/unit_test.hpp>
#include <xtensor/xcomplex.hpp>

#include "common/ParameterSet.h"
#include "steps/ResultStep.h"

using dp3::steps::Clipper;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(clipper)

namespace {
const size_t kNBaselines = 1;
const size_t kNChannels = 4;
const size_t kNCorrelations = 4;
const float kMaxAmplitude = 1.0;
}  // namespace

/// A simple class to test if Clipper configures
/// its internal Predict step correctly
class MockPredict : public Step {
 public:
  MockPredict(int frequency_step) : frequency_step_(frequency_step) {}

  dp3::common::Fields getRequiredFields() const override { return kUvwField; }

  void finish() { getNextStep()->finish(); }
  void show(std::ostream& os) const {}
  dp3::common::Fields getProvidedFields() const override { return kDataField; }

  bool process(std::unique_ptr<dp3::base::DPBuffer> buffer) override {
    const dp3::base::DPBuffer::UvwType& uvw_buffer = buffer->GetUvw();
    const std::array<std::size_t, 2> expected_uvw_shape = {kNBaselines, 3};
    BOOST_TEST(uvw_buffer.shape() == expected_uvw_shape,
               boost::test_tools::per_element());

    dp3::base::DPBuffer::DataType& data = buffer->GetData();
    std::array<size_t, 3> data_shape{kNBaselines, kNChannels, kNCorrelations};
    data.resize(data_shape);
    size_t size = buffer->GetData().size();
    for (size_t i = 0; i < size; i++) {
      // This range results in visibilities with real part between
      // 0 and 2*kMaxAmplitude. Clipper should flag the second half
      // of the visibilities in the test.
      data[i] = std::complex<float>(2 * kMaxAmplitude * frequency_step_ * i /
                                        static_cast<float>(size - 1),
                                    0.0);
    }

    getNextStep()->process(std::move(buffer));
    return false;
  }

  void updateInfo(const dp3::base::DPInfo& info_in) override {
    Step::updateInfo(info_in);
    BOOST_TEST(kNChannels == info().origNChan());
    BOOST_TEST(kNChannels / frequency_step_ == info_in.nchan());
  }

 private:
  const int frequency_step_;
};

void TestClipper(size_t time_step, size_t frequency_step) {
  dp3::common::ParameterSet parset;
  parset.add("timestep", std::to_string(time_step));
  parset.add("freqstep", std::to_string(frequency_step));
  parset.add("amplmax", std::to_string(kMaxAmplitude));
  parset.add("sourcedb", "tNDPPP-generic.MS/sky");

  dp3::base::DPInfo info(kNCorrelations, kNChannels);
  info.setTimes(0.0, 2.0, 1.0);
  info.setChannels(std::vector<double>(kNChannels, 42.0e6),
                   std::vector<double>(kNChannels, 1.0e6));

  auto buffer = std::make_unique<dp3::base::DPBuffer>();
  const std::array<size_t, 2> uvw_shape = {kNBaselines, 3};
  buffer->GetUvw() =
      xt::arange<double>(0.0, kNBaselines * kNChannels * kNCorrelations, 1)
          .reshape(uvw_shape);

  std::array<size_t, 3> data_shape = {kNBaselines, kNChannels, kNCorrelations};
  buffer->GetFlags().resize(data_shape);
  buffer->GetFlags().fill(false);

  auto clipper_step = std::make_shared<Clipper>(parset, "");
  auto result_step = std::make_shared<dp3::steps::ResultStep>();
  auto mock_predict = std::make_shared<MockPredict>(frequency_step);
  clipper_step->setNextStep(result_step);
  clipper_step->SetPredict(mock_predict);
  clipper_step->setInfo(info);

  const xt::xtensor<bool, 3> expected_flags{{{false, false, false, false},
                                             {false, false, false, false},
                                             {true, true, true, true},
                                             {true, true, true, true}}};

  for (size_t t = 0; t < info.ntime(); ++t) {
    clipper_step->process(std::make_unique<dp3::base::DPBuffer>(*buffer));

    BOOST_TEST(result_step->get().GetFlags() == expected_flags,
               boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(test_basic_clipper) {
  TestClipper(1, 1);
  TestClipper(1, 2);
  TestClipper(2, 1);
  TestClipper(2, 2);
}

BOOST_AUTO_TEST_SUITE_END()
