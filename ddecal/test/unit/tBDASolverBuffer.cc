// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/BDASolverBuffer.h"

#include "../../../base/BDABuffer.h"
#include "../../../base/test/unit/tBDABuffer.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

using dp3::base::BDABuffer;
using dp3::base::BDASolverBuffer;

namespace {
const size_t kNBaselines = 3;
const size_t kNChannels = 6;
const size_t kNCorrelations = 4;
const double kStartTime = 42.0;
const size_t kInterval = 2.0;
const double kTolerance = 1e-6;

/**
 * Create successive BDA buffers.
 * The resulting buffers hold values for all fields, including flags, which
 * allows testing if AssignAndWeight discards unnecessary fields.
 * @param data_value Visibility value for all buffer entries.
 * @param weight_value Weight value for all buffer entries.
 * @return A series of successive BDA Buffers.
 */
std::vector<std::unique_ptr<BDABuffer>> CreateBDABuffers(
    std::complex<float> data_value, float weight_value) {
  std::vector<std::unique_ptr<BDABuffer>> buffers;

  const std::vector<std::complex<float>> data(kNChannels * kNCorrelations,
                                              data_value);
  const std::vector<float> weights(kNChannels * kNCorrelations, weight_value);
  const bool flags[kNChannels * kNCorrelations]{false};

  double time = kStartTime;
  for (std::size_t i = 0; i < 3; ++i) {
    auto buffer = boost::make_unique<BDABuffer>(kNBaselines * kNChannels *
                                                kNCorrelations);
    // Add a non-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 0, kNChannels, kNCorrelations,
                   data.data(), flags, weights.data(), flags, nullptr);
    buffer->AddRow(time + kInterval, kInterval, kInterval, 0, kNChannels,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    // Add a time-averaged baseline.
    buffer->AddRow(time, kInterval * 2.0, kInterval * 2.0, 1, kNChannels,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    // Add a frequency-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    buffer->AddRow(time + kInterval, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    buffers.emplace_back(std::move(buffer));
    time += kInterval * 2.0;
  }

  return buffers;
}

void CheckComplex(const std::complex<float>& left,
                  const std::complex<float>& right) {
  BOOST_CHECK_CLOSE(left.real(), right.real(), kTolerance);
  BOOST_CHECK_CLOSE(left.imag(), right.imag(), kTolerance);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bdasolverbuffer)

BOOST_AUTO_TEST_CASE(constructor) {
  const BDASolverBuffer buffer;
  BOOST_CHECK(buffer.GetDataRows().empty());
}

BOOST_AUTO_TEST_CASE(set_directions) {
  const size_t kNDirections = 42;
  BDASolverBuffer buffer;
  buffer.SetDirections(kNDirections);
  for (size_t dir = 0; dir < kNDirections; ++dir) {
    BOOST_CHECK(buffer.GetModelDataRows(dir).empty());
  }
}

BOOST_AUTO_TEST_CASE(append_and_weight) {
  const std::complex<float> kDataValue{1.0, 2.0};
  const std::complex<float> kModelValue1 = {3.0, 4.0};
  const std::complex<float> kModelValue2 = {4.0, 5.0};
  const float kWeight = 0.8;
  const float kWeightSqrt = std::sqrt(kWeight);
  // AssignAndWeight should ignore the model weight values.
  const float kModelWeight = 0.0;
  const std::complex<float> kWeightedDataValue = kDataValue * kWeightSqrt;
  const std::complex<float> kWeightedModelValue1 = kModelValue1 * kWeightSqrt;
  const std::complex<float> kWeightedModelValue2 = kModelValue2 * kWeightSqrt;

  std::vector<std::unique_ptr<BDABuffer>> data =
      CreateBDABuffers(kDataValue, kWeight);
  std::vector<std::unique_ptr<BDABuffer>> model_data1 =
      CreateBDABuffers(kModelValue1, kModelWeight);
  std::vector<std::unique_ptr<BDABuffer>> model_data2 =
      CreateBDABuffers(kModelValue2, kModelWeight);

  BDASolverBuffer buffer;
  buffer.SetDirections(2);
  buffer.SetInterval(0.0, 1.0);
  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<std::unique_ptr<BDABuffer>> model_buffers;
    model_buffers.push_back(std::move(model_data1[i]));
    model_buffers.push_back(std::move(model_data2[i]));
    buffer.AppendAndWeight(*data[i], std::move(model_buffers));
  }

  BOOST_REQUIRE(buffer.GetData().size() == data.size());
  BOOST_REQUIRE(buffer.GetModelData(0).size() == data.size());
  BOOST_REQUIRE(buffer.GetModelData(1).size() == data.size());

  for (size_t b = 0; b < data.size(); ++b) {
    const BDABuffer& data_buffer = buffer.GetData()[b];
    const BDABuffer& model_buffer1 = *buffer.GetModelData(0)[b];
    const BDABuffer& model_buffer2 = *buffer.GetModelData(1)[b];

    // AssignAndWeight should only store the visibilities in the buffer.
    BOOST_CHECK(data_buffer.GetData());
    BOOST_CHECK(!data_buffer.GetFlags());
    BOOST_CHECK(!data_buffer.GetWeights());
    BOOST_CHECK(!data_buffer.GetFullResFlags());

    dp3::base::test::CheckBDARowMetaData(*data[b], data_buffer);
    dp3::base::test::CheckBDARowMetaData(*data[b], model_buffer1);
    dp3::base::test::CheckBDARowMetaData(*data[b], model_buffer2);

    const size_t n_elements = data[b]->GetNumberOfElements();
    BOOST_REQUIRE(data_buffer.GetNumberOfElements() == n_elements);
    BOOST_REQUIRE(model_buffer1.GetNumberOfElements() == n_elements);
    BOOST_REQUIRE(model_buffer2.GetNumberOfElements() == n_elements);
    for (size_t i = 0; i < n_elements; ++i) {
      CheckComplex(data_buffer.GetData()[i], kWeightedDataValue);
      CheckComplex(model_buffer1.GetData()[i], kWeightedModelValue1);
      CheckComplex(model_buffer2.GetData()[i], kWeightedModelValue2);
    }
  }

  buffer.Clear();
  BOOST_CHECK(buffer.GetData().empty());
}

BOOST_AUTO_TEST_SUITE_END()
