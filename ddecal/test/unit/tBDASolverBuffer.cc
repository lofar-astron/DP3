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
const double kFirstTime = 42.0;
const size_t kInterval = 2.0;
const size_t kRowsPerBuffer = 5;
const double kTolerance = 1e-6;
const std::complex<float> kDataIncrement{1.0, 1.0};  // Value increment per row.
const std::complex<float> kFirstDataValue{1.0, 2.0};
const std::complex<float> kModelDataDiff{100.0, 100.0};

/**
 * Create successive BDA buffers.
 * The resulting buffers hold values for all fields, including flags, which
 * allows testing if AssignAndWeight discards unnecessary fields.
 * @param data_value Visibility value for the first row. For each next row, the
 *        value is increased with { 1.0, 1.0 }.
 * @param weight_value Weight value for all buffer entries.
 * @return A series of successive BDA Buffers.
 */
std::vector<std::unique_ptr<BDABuffer>> CreateBDABuffers(
    std::complex<float> data_value, float weight_value) {
  std::vector<std::unique_ptr<BDABuffer>> buffers;

  std::vector<std::complex<float>> data(kNChannels * kNCorrelations,
                                        data_value);
  const std::vector<float> weights(kNChannels * kNCorrelations, weight_value);
  const bool flags[kNChannels * kNCorrelations]{false};

  double time = kFirstTime;
  for (std::size_t i = 0; i < 3; ++i) {
    auto buffer = boost::make_unique<BDABuffer>(kRowsPerBuffer * kNChannels *
                                                kNCorrelations);
    // Add a non-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 0, kNChannels, kNCorrelations,
                   data.data(), flags, weights.data(), flags, nullptr);
    for (std::complex<float>& d : data) d += kDataIncrement;
    buffer->AddRow(time + kInterval, kInterval, kInterval, 0, kNChannels,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    for (std::complex<float>& d : data) d += kDataIncrement;

    // Add a time-averaged baseline.
    buffer->AddRow(time, kInterval * 2.0, kInterval * 2.0, 1, kNChannels,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    for (std::complex<float>& d : data) d += kDataIncrement;

    // Add a frequency-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    for (std::complex<float>& d : data) d += kDataIncrement;

    buffer->AddRow(time + kInterval, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data(), flags,
                   nullptr);
    for (std::complex<float>& d : data) d += kDataIncrement;

    BOOST_REQUIRE_EQUAL(buffer->GetRows().size(), kRowsPerBuffer);
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

/**
 * Compares a data row and the first model row against a reference row.
 */
void CheckRow(const BDASolverBuffer& solver_buffer, size_t row_index,
              const BDABuffer::Row& reference_row) {
  // Get the row from the solver buffer.
  const std::vector<const BDABuffer::Row*>& data_rows =
      solver_buffer.GetDataRows();
  const std::vector<const BDABuffer::Row*>& model_rows =
      solver_buffer.GetModelDataRows(0);
  BOOST_REQUIRE_LT(row_index, data_rows.size());
  BOOST_REQUIRE_EQUAL(data_rows.size(), model_rows.size());
  BOOST_REQUIRE(data_rows[row_index] && model_rows[row_index]);
  const BDABuffer::Row& data_row = *data_rows[row_index];
  const BDABuffer::Row& model_row = *model_rows[row_index];

  // Compare the solver buffer row against the original row.
  dp3::base::test::CheckBDARowMetaData(reference_row, data_row);
  dp3::base::test::CheckBDARowMetaData(reference_row, model_row);
  // Checking the value for the first correlation only should suffice.
  CheckComplex(*data_row.data, *reference_row.data);
  CheckComplex(*model_row.data, *reference_row.data + kModelDataDiff);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bdasolverbuffer)

BOOST_AUTO_TEST_CASE(constructor) {
  const BDASolverBuffer buffer(1, 0.0, 1.0);
  BOOST_CHECK_EQUAL(buffer.BufferCount(), 0u);
  BOOST_CHECK(buffer.GetDataRows().empty());
  BOOST_CHECK(buffer.GetModelDataRows(0).empty());
}

BOOST_AUTO_TEST_CASE(append_and_weight) {
  const float kWeight = 0.8;
  const float kWeightSqrt = std::sqrt(kWeight);
  // AssignAndWeight should ignore the model weight values.
  const float kModelWeight = 0.0;

  std::vector<std::unique_ptr<BDABuffer>> data =
      CreateBDABuffers(kFirstDataValue, kWeight);
  std::vector<std::unique_ptr<BDABuffer>> model_data1 =
      CreateBDABuffers(kFirstDataValue + 1.0f * kModelDataDiff, kModelWeight);
  std::vector<std::unique_ptr<BDABuffer>> model_data2 =
      CreateBDABuffers(kFirstDataValue + 2.0f * kModelDataDiff, kModelWeight);

  // All bda rows are in the first interval.
  BDASolverBuffer buffer(2, 0.0, 100.0);
  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<std::unique_ptr<BDABuffer>> model_buffers;
    model_buffers.push_back(std::move(model_data1[i]));
    model_buffers.push_back(std::move(model_data2[i]));
    buffer.AppendAndWeight(*data[i], std::move(model_buffers));
  }

  const size_t n_rows = data.size() * kRowsPerBuffer;
  BOOST_REQUIRE_EQUAL(buffer.GetDataRows().size(), n_rows);
  BOOST_REQUIRE_EQUAL(buffer.GetModelDataRows(0).size(), n_rows);
  BOOST_REQUIRE_EQUAL(buffer.GetModelDataRows(1).size(), n_rows);

  for (size_t row = 0; row < n_rows; ++row) {
    const BDABuffer::Row* data_row = buffer.GetDataRows()[row];
    const BDABuffer::Row* model_row1 = buffer.GetModelDataRows(0)[row];
    const BDABuffer::Row* model_row2 = buffer.GetModelDataRows(1)[row];

    BOOST_REQUIRE(data_row);
    BOOST_REQUIRE(model_row1);
    BOOST_REQUIRE(model_row2);

    // AssignAndWeight should only store the visibilities in the buffer.
    BOOST_CHECK(data_row->data);
    BOOST_CHECK(!data_row->flags);
    BOOST_CHECK(!data_row->weights);
    BOOST_CHECK(!data_row->full_res_flags);

    const size_t data_buffer_index = row / kRowsPerBuffer;
    const size_t data_row_index = row % kRowsPerBuffer;
    const BDABuffer::Row& original_row =
        data[data_buffer_index]->GetRows()[data_row_index];

    dp3::base::test::CheckBDARowMetaData(original_row, *data_row);
    dp3::base::test::CheckBDARowMetaData(original_row, *model_row1);
    dp3::base::test::CheckBDARowMetaData(original_row, *model_row2);

    for (size_t i = 0; i < original_row.GetDataSize(); ++i) {
      const std::complex<float> inc = float(row) * kDataIncrement;
      const std::complex<float> weighted_data_value =
          (kFirstDataValue + inc) * kWeightSqrt;
      const std::complex<float> weighted_model_value1 =
          (kFirstDataValue + 1.0f * kModelDataDiff + inc) * kWeightSqrt;
      const std::complex<float> weighted_model_value2 =
          (kFirstDataValue + 2.0f * kModelDataDiff + inc) * kWeightSqrt;
      CheckComplex(data_row->data[i], weighted_data_value);
      CheckComplex(model_row1->data[i], weighted_model_value1);
      CheckComplex(model_row2->data[i], weighted_model_value2);
    }
  }

  buffer.Clear();
  BOOST_CHECK_EQUAL(buffer.GetDataRows().size(), 0u);
}

BOOST_AUTO_TEST_CASE(one_interval_per_buffer) {
  const float kWeight = 1.0;

  std::vector<std::unique_ptr<BDABuffer>> data =
      CreateBDABuffers(kFirstDataValue, kWeight);

  BDASolverBuffer buffer(0, kFirstTime - kInterval / 2, kInterval * 2.0);
  for (size_t i = 0; i < data.size(); ++i) {
    buffer.AppendAndWeight(*data[i], {});
    BOOST_CHECK_EQUAL(buffer.BufferCount(), i + 1);
  }

  for (size_t data_index = 0; data_index < data.size(); ++data_index) {
    const bool complete = data_index + 1 < data.size();
    BOOST_CHECK_EQUAL(complete, buffer.IntervalIsComplete());

    BOOST_REQUIRE_EQUAL(buffer.GetDataRows().size(), kRowsPerBuffer);
    for (size_t row_index = 0; row_index < kRowsPerBuffer; ++row_index) {
      const BDABuffer::Row* row = buffer.GetDataRows()[row_index];
      BOOST_REQUIRE(row && row->data);

      const BDABuffer::Row& original_row =
          data[data_index]->GetRows()[row_index];

      dp3::base::test::CheckBDARowMetaData(original_row, *row);

      for (size_t i = 0; i < original_row.GetDataSize(); ++i) {
        const size_t all_row_index = data_index * kRowsPerBuffer + row_index;
        const std::complex<float> weighted_data_value =
            kFirstDataValue + float(all_row_index) * kDataIncrement;
        CheckComplex(row->data[i], weighted_data_value);
      }
    }
    buffer.AdvanceInterval();
  }

  BOOST_CHECK_EQUAL(buffer.GetDataRows().size(), 0u);
}

BOOST_AUTO_TEST_CASE(buffer_with_multiple_intervals) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;

  // Create a BDA buffer with successive rows and a single value per row.
  BDABuffer data_buffer(kNRows * kNCorrelations);
  auto model_buffer = boost::make_unique<BDABuffer>(kNRows * kNCorrelations);

  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0};
  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffer.AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                       data_value.data(), nullptr, kWeight.data());
    model_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                         model_value.data(), nullptr, kWeight.data());
    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  BDASolverBuffer solver_buffer(1, kFirstTime - kInterval / 2,
                                kSolutionInterval);
  std::vector<std::unique_ptr<BDABuffer>> model_buffers;
  model_buffers.push_back(std::move(model_buffer));
  solver_buffer.AppendAndWeight(data_buffer, std::move(model_buffers));
  BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 1u);

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, row_index, data_buffer.GetRows()[row]);

    ++row_index;
    if (row_index == kSolutionIntervalFactor) {
      row_index = 0;
      solver_buffer.AdvanceInterval();
      const bool complete = row < kNRows - kSolutionIntervalFactor;
      BOOST_CHECK_EQUAL(solver_buffer.IntervalIsComplete(), complete);
      BOOST_CHECK_EQUAL(solver_buffer.GetDataRows().size(),
                        complete ? kSolutionIntervalFactor
                                 : (kNRows % kSolutionIntervalFactor));
      BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 1u);
    }
  }
}

BOOST_AUTO_TEST_CASE(multiple_buffers_per_interval) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;

  // Create successive BDABuffers with a single row.
  std::vector<std::unique_ptr<BDABuffer>> data_buffers;
  std::vector<std::unique_ptr<BDABuffer>> model_buffers;
  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0};

  BDASolverBuffer solver_buffer(1, kFirstTime - kInterval / 2,
                                kSolutionInterval);

  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffers.push_back(boost::make_unique<BDABuffer>(kNCorrelations));
    auto model_buffer = boost::make_unique<BDABuffer>(kNCorrelations);
    BOOST_REQUIRE(data_buffers.back()->AddRow(time, kInterval, kInterval, 0, 1,
                                              kNCorrelations, data_value.data(),
                                              nullptr, kWeight.data()));
    BOOST_REQUIRE(model_buffer->AddRow(time, kInterval, kInterval, 0, 1,
                                       kNCorrelations, model_value.data(),
                                       nullptr, kWeight.data()));

    std::vector<std::unique_ptr<BDABuffer>> model_buffers;
    model_buffers.push_back(std::move(model_buffer));
    solver_buffer.AppendAndWeight(*data_buffers.back(),
                                  std::move(model_buffers));
    BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), row + 1);

    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, row_index, data_buffers[row]->GetRows().front());

    ++row_index;
    if (row_index == kSolutionIntervalFactor) {
      row_index = 0;
      solver_buffer.AdvanceInterval();
      const bool complete = row < kNRows - kSolutionIntervalFactor;
      BOOST_CHECK_EQUAL(solver_buffer.IntervalIsComplete(), complete);
      BOOST_CHECK_EQUAL(solver_buffer.GetDataRows().size(),
                        complete ? kSolutionIntervalFactor
                                 : (kNRows % kSolutionIntervalFactor));
      BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), kNRows - row - 1);
    }
  }

  // Test advancing beyond the last interval.
  solver_buffer.AdvanceInterval();
  BOOST_CHECK(!solver_buffer.IntervalIsComplete());
  BOOST_CHECK_EQUAL(solver_buffer.GetDataRows().size(), 0u);
  BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()
