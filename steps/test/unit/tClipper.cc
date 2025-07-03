// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/Clipper.h"

#include <algorithm>
#include <iostream>
#include <filesystem>
#include <fstream>

#include <boost/test/unit_test.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xview.hpp>

#include "common/ParameterSet.h"
#include "steps/ResultStep.h"
#include "tPredict.h"

using dp3::steps::Clipper;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(clipper)

namespace {
const size_t kNBaselines = 1;
const size_t kNChannels = 4;
const size_t kNCorrelations = 4;
const float kMaxAmplitude = 1.0;
const size_t kFlagDistance = std::min(kNChannels, kNCorrelations);
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
      // This parametrization results in a visibility amplitude which (for each
      // baseline) increases linearly with the distance from the first (channel,
      // correlation) bin. Visibilities exceed kMaxAmplitude in a sawtooth-like
      // structure. E.g. for kNBaselines = 1, kNChannels = kNCorrelations = 40,
      // kMaxAmplitude = 4, one has data =
      //        {{{(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)},
      //          {(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)},
      //          {(2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0)},
      //          {(3.0, 0.0), (4.0, 0.0), (5.0, 0.0), (6.0, 0.0)}}}
      // The visibilities in the lower right exceed kMaxAmplitude.
      data[i] = std::complex<float>(
          kMaxAmplitude *
              (i / kNCorrelations * frequency_step_ + i % kNCorrelations) /
              static_cast<float>(kFlagDistance),
          0.0);
    }

    getNextStep()->process(std::move(buffer));
    return false;
  }

  void updateInfo(const dp3::base::DPInfo& info_in) override {
    Step::updateInfo(info_in);
    BOOST_TEST(kNChannels == getInfoOut().origNChan());
    BOOST_TEST(kNChannels / frequency_step_ == info_in.nchan());
  }

 private:
  const int frequency_step_;
};

void TestClipper(size_t time_step, size_t frequency_step,
                 bool all_correlations) {
  dp3::common::ParameterSet parset;
  parset.add("timestep", std::to_string(time_step));
  parset.add("freqstep", std::to_string(frequency_step));
  parset.add("amplmax", std::to_string(kMaxAmplitude));
  parset.add("flagallcorrelations", std::to_string(all_correlations));
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);

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

  auto expected_flags =
      xt::xarray<bool>::from_shape({kNBaselines, kNChannels, kNCorrelations});
  expected_flags.fill(false);
  const size_t size = expected_flags.size();
  // The expected flags assume the sawtooth structure described above,
  // except the "predicted" values are repeated `frequency_step` per channel.
  for (size_t i = 0; i < size; i++) {
    const size_t distance_from_start = i / kNCorrelations + i % kNCorrelations;
    const size_t n_channels_old_predict = (i / kNCorrelations) % frequency_step;
    if (distance_from_start > kFlagDistance) expected_flags[i] = true;
    if (n_channels_old_predict)
      expected_flags[i] =
          expected_flags[i - n_channels_old_predict * kNCorrelations];
  }
  xt::xarray<int> contains_true = xt::sum(expected_flags, -1);
  if (all_correlations) {
    auto idx = xt::from_indices(xt::argwhere(contains_true));
    auto selection = xt::view(idx, xt::all(), 1);
    auto view =
        xt::view(expected_flags, xt::all(), xt::keep(selection), xt::all());
    view = true;
  }

  for (size_t t = 0; t < info.ntime(); ++t) {
    clipper_step->process(std::make_unique<dp3::base::DPBuffer>(*buffer));

    BOOST_TEST(result_step->get().GetFlags() == expected_flags,
               boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(test_basic_clipper) {
  for (size_t time_step : {1, 2}) {
    for (size_t frequency_step : {1, 2}) {
      for (bool all_correlations : {true, false}) {
        TestClipper(time_step, frequency_step, all_correlations);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
