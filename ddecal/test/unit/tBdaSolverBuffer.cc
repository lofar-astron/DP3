// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/BdaSolverBuffer.h"

#include <dp3/base/BdaBuffer.h>

#include <boost/test/unit_test.hpp>

#include <aocommon/matrix2x2.h>

#include <array>

using dp3::base::BdaBuffer;
using dp3::ddecal::BdaSolverBuffer;

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
 * allows testing if AppendAndWeight discards unnecessary fields.
 * @param data_value Visibility value for the first row. For each next row, the
 *        value is increased with { 1.0, 1.0 }.
 * @param weight_value Weight value for all buffer entries.
 * @return A series of successive BDA Buffers.
 */
std::vector<std::unique_ptr<BdaBuffer>> CreateBdaBuffers(
    std::complex<float> data_value, float weight_value) {
  std::vector<std::unique_ptr<BdaBuffer>> buffers;

  std::vector<std::complex<float>> data(kNChannels * kNCorrelations,
                                        data_value);
  const std::vector<float> weights(kNChannels * kNCorrelations, weight_value);
  const bool flags[kNChannels * kNCorrelations]{false};

  double time = kFirstTime;
  for (std::size_t i = 0; i < 3; ++i) {
    auto buffer = std::make_unique<BdaBuffer>(kRowsPerBuffer * kNChannels *
                                              kNCorrelations);
    // Add a non-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 0, kNChannels, kNCorrelations,
                   data.data(), flags, weights.data());
    for (std::complex<float>& d : data) d += kDataIncrement;
    buffer->AddRow(time + kInterval, kInterval, kInterval, 0, kNChannels,
                   kNCorrelations, data.data(), flags, weights.data());
    for (std::complex<float>& d : data) d += kDataIncrement;

    // Add a time-averaged baseline.
    buffer->AddRow(time + kInterval / 2, kInterval * 2.0, kInterval * 2.0, 1,
                   kNChannels, kNCorrelations, data.data(), flags,
                   weights.data());
    for (std::complex<float>& d : data) d += kDataIncrement;

    // Add a frequency-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data());
    for (std::complex<float>& d : data) d += kDataIncrement;

    buffer->AddRow(time + kInterval, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations, data.data(), flags, weights.data());
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
 * Compares the meta data of a BdaBuffer::Row and an IntervalRow.
 */
void CheckRowMetaData(const BdaBuffer::Row& bda_row,
                      const BdaSolverBuffer::IntervalRow& solver_row) {
  BOOST_TEST(bda_row.time == solver_row.time);
  BOOST_TEST(bda_row.baseline_nr == solver_row.baseline_nr);
  BOOST_TEST(bda_row.n_channels == solver_row.n_channels);
  BOOST_TEST(bda_row.n_correlations == solver_row.n_correlations);
}

/**
 * Compares a data row and the first model row against a reference row.
 */
void CheckRow(const BdaSolverBuffer& solver_buffer, size_t row_index,
              const BdaBuffer::Row& reference_row,
              const std::complex<float>* reference_row_data) {
  // Get the row from the solver buffer.
  const std::vector<BdaSolverBuffer::IntervalRow>& rows =
      solver_buffer.GetIntervalRows();
  BOOST_REQUIRE_LT(row_index, rows.size());
  const BdaSolverBuffer::IntervalRow& row = rows[row_index];
  BOOST_REQUIRE(row.weighted_data);
  BOOST_REQUIRE(!row.model_data.empty());
  BOOST_REQUIRE(row.model_data.front());

  // Compare the solver buffer row against the original row.
  CheckRowMetaData(reference_row, row);
  // Checking the value for the first correlation only should suffice.
  CheckComplex(*row.weighted_data, *reference_row_data);
  CheckComplex(*row.model_data.front(), *reference_row_data + kModelDataDiff);
}

/**
 * Utility function for calling AppendAndWeight with a single model buffer.
 */
void DoAppendAndWeight(BdaSolverBuffer& solver_buffer,
                       std::unique_ptr<BdaBuffer> data_buffer,
                       std::unique_ptr<BdaBuffer> model_buffer) {
  std::vector<std::unique_ptr<BdaBuffer>> model_buffers;
  model_buffers.push_back(std::move(model_buffer));
  solver_buffer.AppendAndWeight(std::move(data_buffer),
                                std::move(model_buffers));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bdasolverbuffer)

BOOST_AUTO_TEST_CASE(constructor) {
  BdaSolverBuffer buffer(1, 0.0, 1.0, 1);
  BOOST_CHECK_EQUAL(buffer.BufferCount(), 0u);
  BOOST_CHECK(buffer.GetIntervalRows().empty());
  BOOST_CHECK(buffer.GetDone().empty());
}

BOOST_AUTO_TEST_CASE(append_and_weight) {
  const float kWeight = 0.8;
  const float kWeightSqrt = std::sqrt(kWeight);
  // AppendAndWeight should ignore the model weight values.
  const float kModelWeight = 0.0;

  std::vector<std::unique_ptr<BdaBuffer>> data =
      CreateBdaBuffers(kFirstDataValue, kWeight);
  std::vector<BdaBuffer*> data_ptr;
  std::vector<std::unique_ptr<BdaBuffer>> model_data1 =
      CreateBdaBuffers(kFirstDataValue + 1.0f * kModelDataDiff, kModelWeight);
  std::vector<std::unique_ptr<BdaBuffer>> model_data2 =
      CreateBdaBuffers(kFirstDataValue + 2.0f * kModelDataDiff, kModelWeight);

  // All bda rows are in the first interval.
  BdaSolverBuffer buffer(2, 0.0, 100.0, kNBaselines);
  for (size_t i = 0; i < data.size(); ++i) {
    data_ptr.push_back(data[i].get());
    std::vector<std::unique_ptr<BdaBuffer>> model_buffers;
    model_buffers.push_back(std::move(model_data1[i]));
    model_buffers.push_back(std::move(model_data2[i]));
    buffer.AppendAndWeight(std::move(data[i]), std::move(model_buffers));
  }

  const size_t n_rows = data.size() * kRowsPerBuffer;
  BOOST_REQUIRE_EQUAL(buffer.GetIntervalRows().size(), n_rows);

  for (size_t row_index = 0; row_index < n_rows; ++row_index) {
    const BdaSolverBuffer::IntervalRow& solver_row =
        buffer.GetIntervalRows()[row_index];
    BOOST_REQUIRE(solver_row.weighted_data);
    BOOST_REQUIRE_EQUAL(solver_row.model_data.size(), 2);
    BOOST_REQUIRE(solver_row.model_data[0]);
    BOOST_REQUIRE(solver_row.model_data[1]);

    const size_t data_buffer_index = row_index / kRowsPerBuffer;
    const size_t data_row_index = row_index % kRowsPerBuffer;
    const BdaBuffer::Row& original_row =
        data_ptr[data_buffer_index]->GetRows()[data_row_index];

    CheckRowMetaData(original_row, solver_row);

    for (size_t i = 0; i < original_row.GetDataSize(); ++i) {
      const std::complex<float> inc = float(row_index) * kDataIncrement;
      const std::complex<float> weighted_data_value =
          (kFirstDataValue + inc) * kWeightSqrt;
      const std::complex<float> weighted_model_value1 =
          (kFirstDataValue + 1.0f * kModelDataDiff + inc) * kWeightSqrt;
      const std::complex<float> weighted_model_value2 =
          (kFirstDataValue + 2.0f * kModelDataDiff + inc) * kWeightSqrt;
      CheckComplex(solver_row.weighted_data[i], weighted_data_value);
      CheckComplex(solver_row.model_data[0][i], weighted_model_value1);
      CheckComplex(solver_row.model_data[1][i], weighted_model_value2);
    }
  }

  buffer.Clear();
  BOOST_TEST(buffer.GetIntervalRows().empty());
}

BOOST_AUTO_TEST_CASE(append_and_weight_nullptr_flags) {
  const std::vector<std::complex<float>> data(kNChannels * kNCorrelations, 1.0);
  const std::vector<float> weights(kNChannels * kNCorrelations, 1.0);

  BdaBuffer::Fields fields(false);
  fields.data = true;
  fields.weights = true;
  auto bda_buffer =
      std::make_unique<BdaBuffer>(kNChannels * kNCorrelations, fields);
  bda_buffer->AddRow(kFirstTime, kInterval, kInterval, 0, kNChannels,
                     kNCorrelations, data.data(), nullptr, weights.data());

  auto model_data =
      std::make_unique<BdaBuffer>(kRowsPerBuffer * kNChannels * kNCorrelations);
  model_data->AddRow(kFirstTime, kInterval, kInterval, 0, kNChannels,
                     kNCorrelations, data.data(), nullptr, weights.data());

  BdaSolverBuffer buffer(1, 0.0, 100.0, kNBaselines);
  std::vector<std::unique_ptr<BdaBuffer>> model_buffers;
  model_buffers.push_back(std::move(model_data));
  BOOST_CHECK_NO_THROW(
      buffer.AppendAndWeight(std::move(bda_buffer), std::move(model_buffers)));
}

BOOST_AUTO_TEST_CASE(one_interval_per_buffer) {
  const float kWeight = 1.0;

  std::vector<std::unique_ptr<BdaBuffer>> data =
      CreateBdaBuffers(kFirstDataValue, kWeight);
  std::vector<BdaBuffer*> data_ptr;

  BdaSolverBuffer buffer(0, kFirstTime - kInterval / 2, kInterval * 2.0,
                         kNBaselines);
  for (size_t i = 0; i < data.size(); ++i) {
    data_ptr.push_back(data[i].get());
    buffer.AppendAndWeight(std::move(data[i]), {});
    BOOST_CHECK_EQUAL(buffer.BufferCount(), i + 1);
  }

  for (size_t data_index = 0; data_index < data.size(); ++data_index) {
    const bool complete = data_index + 1 < data.size();
    BOOST_CHECK_EQUAL(complete, buffer.IntervalIsComplete());

    BOOST_REQUIRE_EQUAL(buffer.GetIntervalRows().size(), kRowsPerBuffer);
    for (size_t row_index = 0; row_index < kRowsPerBuffer; ++row_index) {
      const BdaSolverBuffer::IntervalRow& solver_row =
          buffer.GetIntervalRows()[row_index];
      BOOST_REQUIRE(solver_row.weighted_data);

      const BdaBuffer::Row& original_row =
          data_ptr[data_index]->GetRows()[row_index];

      CheckRowMetaData(original_row, solver_row);

      for (size_t i = 0; i < original_row.GetDataSize(); ++i) {
        const size_t all_row_index = data_index * kRowsPerBuffer + row_index;
        const std::complex<float> weighted_data_value =
            kFirstDataValue + float(all_row_index) * kDataIncrement;
        CheckComplex(solver_row.weighted_data[i], weighted_data_value);
      }
    }
    buffer.AdvanceInterval();

    const std::vector<std::unique_ptr<BdaBuffer>> done = buffer.GetDone();
    BOOST_REQUIRE_EQUAL(done.size(), 1u);
    BOOST_CHECK_EQUAL(done.front().get(), data_ptr[data_index]);
  }

  BOOST_CHECK(buffer.GetIntervalRows().empty());
}

BOOST_AUTO_TEST_CASE(buffer_with_multiple_intervals) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;

  // Create a BDA buffer with successive rows and a single value per row.
  auto data_buffer = std::make_unique<BdaBuffer>(kNRows * kNCorrelations);
  const BdaBuffer& data_ref = *data_buffer;
  auto model_buffer = std::make_unique<BdaBuffer>(kNRows * kNCorrelations);

  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0, 1.0, 1.0, 1.0};
  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                        data_value.data(), nullptr, kWeight.data());
    model_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                         model_value.data(), nullptr, kWeight.data());
    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  BdaSolverBuffer solver_buffer(1, kFirstTime - kInterval / 2,
                                kSolutionInterval, 1);
  DoAppendAndWeight(solver_buffer, std::move(data_buffer),
                    std::move(model_buffer));
  BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 1u);

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, row_index, data_ref.GetRows()[row],
             data_ref.GetData(row));

    ++row_index;
    if (row_index == kSolutionIntervalFactor) {
      row_index = 0;
      solver_buffer.AdvanceInterval();
      const bool complete = row < kNRows - kSolutionIntervalFactor;
      BOOST_CHECK_EQUAL(solver_buffer.IntervalIsComplete(), complete);
      BOOST_CHECK_EQUAL(solver_buffer.GetIntervalRows().size(),
                        complete ? kSolutionIntervalFactor
                                 : (kNRows % kSolutionIntervalFactor));
      BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 1u);
    }
  }

  // Test that GetDone() returns the data buffer after advancing.
  BOOST_REQUIRE(solver_buffer.GetDone().empty());
  solver_buffer.AdvanceInterval();
  const std::vector<std::unique_ptr<BdaBuffer>> done = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(done.size(), 1u);
  BOOST_CHECK_EQUAL(done.front().get(), &data_ref);
}

BOOST_AUTO_TEST_CASE(multiple_buffers_per_interval) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;

  // Create successive BdaBuffers with a single row.
  std::vector<BdaBuffer*> data_buffers;
  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0, 1.0, 1.0, 1.0};

  BdaSolverBuffer solver_buffer(1, kFirstTime - kInterval / 2,
                                kSolutionInterval, 1);

  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    auto data_buffer = std::make_unique<BdaBuffer>(kNCorrelations);
    data_buffers.push_back(data_buffer.get());
    auto model_buffer = std::make_unique<BdaBuffer>(kNCorrelations);
    BOOST_REQUIRE(data_buffer->AddRow(time, kInterval, kInterval, 0, 1,
                                      kNCorrelations, data_value.data(),
                                      nullptr, kWeight.data()));
    BOOST_REQUIRE(model_buffer->AddRow(time, kInterval, kInterval, 0, 1,
                                       kNCorrelations, model_value.data(),
                                       nullptr, kWeight.data()));

    DoAppendAndWeight(solver_buffer, std::move(data_buffer),
                      std::move(model_buffer));
    BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), row + 1);

    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, row_index, data_buffers[row]->GetRows().front(),
             data_buffers[row]->GetData(0));

    ++row_index;
    if (row_index == kSolutionIntervalFactor) {
      row_index = 0;
      solver_buffer.AdvanceInterval();
      const bool complete = row < kNRows - kSolutionIntervalFactor;
      BOOST_CHECK_EQUAL(solver_buffer.IntervalIsComplete(), complete);
      BOOST_CHECK_EQUAL(solver_buffer.GetIntervalRows().size(),
                        complete ? kSolutionIntervalFactor
                                 : (kNRows % kSolutionIntervalFactor));
      BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), kNRows - row - 1);
    }
  }

  // Test advancing beyond the last interval.
  solver_buffer.AdvanceInterval();
  BOOST_CHECK(!solver_buffer.IntervalIsComplete());
  BOOST_CHECK(solver_buffer.GetIntervalRows().empty());
  BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 0u);

  // Test that GetDone() returns all original data buffers.
  const std::vector<std::unique_ptr<BdaBuffer>> done = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(done.size(), data_buffers.size());
  for (size_t i = 0; i < data_buffers.size(); ++i) {
    BOOST_CHECK_EQUAL(done[i].get(), data_buffers[i]);
  }
}

BOOST_AUTO_TEST_CASE(subtractcorrectedmodel) {
  // Test full jones subtraction on the first row,
  // test diagonal subtraction on the second row,
  // test scalar subtraction on the third row,
  // and test that the fourth row remains intact.
  const size_t kNRows = 4;
  const aocommon::MC2x2F kDataValue({10, 20}, {11, 21}, {12, 22}, {13, 23});
  const aocommon::MC2x2F kModelValue({1, 10}, {2, 20}, {3, 30}, {4, 40});
  const std::array<float, kNCorrelations> kWeight{1, 1, 1, 1};
  const std::vector<std::vector<std::complex<double>>> kSolutions{
      {{1, 2}, {2, 2}, {3, 4}, {4, 4}, {5, 2}, {6, 2}, {7, 1}, {8, 1}}};

  const aocommon::MC2x2F kSolution1(kSolutions.front().data());
  const aocommon::MC2x2F kSolution2(kSolutions.front().data() + kNCorrelations);
  const aocommon::MC2x2F kSolvedFullJones =
      kSolution1.Multiply(kModelValue).MultiplyHerm(kSolution2);
  const aocommon::MC2x2F kFullJonesResult(kDataValue - kSolvedFullJones);

  const aocommon::MC2x2F kDiagonalSolution1(
      std::complex<float>(kSolutions.front()[0]), 0.0f, 0.0f,
      std::complex<float>(kSolutions.front()[1]));
  const aocommon::MC2x2F kDiagonalSolution2(
      std::complex<float>(kSolutions.front()[2]), 0.0f, 0.0f,
      std::complex<float>(kSolutions.front()[3]));
  const aocommon::MC2x2F kSolvedDiagonal =
      kDiagonalSolution1.Multiply(kModelValue).MultiplyHerm(kDiagonalSolution2);
  const aocommon::MC2x2F kDiagonalResult(kDataValue - kSolvedDiagonal);

  const std::complex<float> kScalarFactor =
      std::complex<float>(kSolutions.front()[0]) *
      std::conj(std::complex<float>(kSolutions.front()[1]));
  const std::array<std::complex<float>, kNCorrelations> kScalarResult{
      kDataValue.Get(0) - kScalarFactor * kModelValue.Get(0),
      kDataValue.Get(1) - kScalarFactor * kModelValue.Get(1),
      kDataValue.Get(2) - kScalarFactor * kModelValue.Get(2),
      kDataValue.Get(3) - kScalarFactor * kModelValue.Get(3)};

  const std::vector<double> kStartFreqs{100, 200};
  const std::vector<std::vector<double>> kChanFreqs{{150}};
  const std::vector<int> kAnt1{0};
  const std::vector<int> kAnt2{1};

  // Create BDA input buffer with two intervals.
  auto data_buffer = std::make_unique<BdaBuffer>(kNRows * kNCorrelations);
  auto model_buffer = std::make_unique<BdaBuffer>(kNRows * kNCorrelations);

  std::complex<float> kDataValue_array[4];
  kDataValue.AssignTo(kDataValue_array);
  std::complex<float> kModelValue_array[4];
  kModelValue.AssignTo(kModelValue_array);

  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                        kDataValue_array, nullptr, kWeight.data());
    model_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                         kModelValue_array, nullptr, kWeight.data());
    time += kInterval;
  }

  BdaSolverBuffer solver_buffer(1, kFirstTime - kInterval / 2, kInterval,
                                kNBaselines);
  DoAppendAndWeight(solver_buffer, std::move(data_buffer),
                    std::move(model_buffer));
  solver_buffer.SubtractCorrectedModel(kSolutions, kStartFreqs, 4, kAnt1, kAnt2,
                                       kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.SubtractCorrectedModel(kSolutions, kStartFreqs, 2, kAnt1, kAnt2,
                                       kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.SubtractCorrectedModel(kSolutions, kStartFreqs, 1, kAnt1, kAnt2,
                                       kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.AdvanceInterval();

  std::vector<std::unique_ptr<BdaBuffer>> data_out = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(data_out.size(), 1u);
  BOOST_REQUIRE_EQUAL(data_out.front()->GetRows().size(), kNRows);
  const std::complex<float>* data_fulljones = data_out.front()->GetData(0);
  const std::complex<float>* data_diagonal = data_out.front()->GetData(1);
  const std::complex<float>* data_scalar = data_out.front()->GetData(2);
  const std::complex<float>* data_uncorrected = data_out.front()->GetData(3);
  for (size_t c = 0; c < kNCorrelations; ++c) {
    CheckComplex(data_fulljones[c], kFullJonesResult.Get(c));
    CheckComplex(data_diagonal[c], kDiagonalResult.Get(c));
    CheckComplex(data_scalar[c], kScalarResult[c]);
    BOOST_CHECK_EQUAL(data_uncorrected[c], kDataValue.Get(c));
  }
}

BOOST_AUTO_TEST_CASE(solution_interval_is_complete) {
  const size_t kSolutionIntervalFactor = 4;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;
  double time = kFirstTime;

  // The first solution interval defined is [41 , 49]
  BdaSolverBuffer solver_buffer(1, kFirstTime - kInterval / 2,
                                kSolutionInterval, 2);
  std::vector<std::complex<float>> data(kNChannels * kNCorrelations, 1.0);
  const std::vector<float> weights(kNChannels * kNCorrelations, 1.0);
  const bool flags[kNChannels * kNCorrelations]{false};

  // 1. Add intervals within the first solution interval
  auto buffer_within_first_solint =
      std::make_unique<BdaBuffer>(6 * kNChannels * kNCorrelations);
  // Baseline 0 - 1: [41, 43]
  buffer_within_first_solint->AddRow(time, kInterval, kInterval, 0, kNChannels,
                                     kNCorrelations, data.data(), flags,
                                     weights.data());
  buffer_within_first_solint->AddRow(time, kInterval, kInterval, 1, kNChannels,
                                     kNCorrelations, data.data(), flags,
                                     weights.data());
  // Baseline 0 - 1: [43, 45]
  buffer_within_first_solint->AddRow(time + kInterval, kInterval, kInterval, 0,
                                     kNChannels, kNCorrelations, data.data(),
                                     flags, weights.data());
  buffer_within_first_solint->AddRow(time + kInterval, kInterval, kInterval, 1,
                                     kNChannels, kNCorrelations, data.data(),
                                     flags, weights.data());
  // Baseline 0: [45, 47]
  buffer_within_first_solint->AddRow(time + 2 * kInterval, kInterval, kInterval,
                                     0, kNChannels, kNCorrelations, data.data(),
                                     flags, weights.data());
  // Baseline 0: [47, 49]
  buffer_within_first_solint->AddRow(time + 3 * kInterval, kInterval, kInterval,
                                     0, kNChannels, kNCorrelations, data.data(),
                                     flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_within_first_solint), {});
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 2. Add intervals on the edge of the solution interval
  auto buffer_edge_solint =
      std::make_unique<BdaBuffer>(2 * kNChannels * kNCorrelations);
  //  Baseline 0: [49, 51]
  buffer_edge_solint->AddRow(time + 4 * kInterval, kInterval, kInterval, 0,
                             kNChannels, kNCorrelations, data.data(), flags,
                             weights.data());
  //  Baseline 1: [45, 51]
  buffer_edge_solint->AddRow(time + 2.5 * kInterval, 3 * kInterval,
                             3 * kInterval, 0, kNChannels, kNCorrelations,
                             data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_edge_solint), {});
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 2. Add interval completely contained in the next solution interval for
  // baseline 0
  auto buffer_within_second_solint_bl_0 =
      std::make_unique<BdaBuffer>(1 * kNChannels * kNCorrelations);
  // Baseline 0: [51, 53]
  buffer_within_second_solint_bl_0->AddRow(
      time + 5 * kInterval, kInterval, kInterval, 0, kNChannels, kNCorrelations,
      data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_within_second_solint_bl_0),
                                {});
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 3. Add interval completely contained in the next solution interval for
  // baseline 1
  auto buffer_within_second_solint_bl_1 =
      std::make_unique<BdaBuffer>(1 * kNChannels * kNCorrelations);
  // Baseline 1: [51, 53]
  buffer_within_second_solint_bl_1->AddRow(
      time + 5 * kInterval, kInterval, kInterval, 1, kNChannels, kNCorrelations,
      data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_within_second_solint_bl_1),
                                {});
  BOOST_CHECK_EQUAL(true, solver_buffer.IntervalIsComplete());
}

BOOST_AUTO_TEST_SUITE_END()
