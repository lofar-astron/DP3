// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/BdaSolverBuffer.h"

#include <array>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <aocommon/matrix2x2.h>

#include <dp3/base/BdaBuffer.h>

#include "../../../common/test/unit/tCommon.h"

using dp3::base::BdaBuffer;
using dp3::common::Fields;
using dp3::common::test::kTrueFalseRange;
using dp3::ddecal::BdaSolverBuffer;

namespace {
const size_t kNBaselines = 3;
const size_t kNChannels = 6;
const size_t kNCorrelations = 4;
const size_t kNElementsPerRow = kNChannels * kNCorrelations;
const double kFirstTime = 42.0;
const double kInterval = 2.0;
const size_t kRowsPerBuffer = 5;
const double kTolerance = 1.0e-6;
const std::complex<float> kDataIncrement{1.0, 1.0};  // Value increment per row.
const std::complex<float> kFirstDataValue{1.0, 2.0};
const std::complex<float> kModelDataDiff{100.0, 100.0};
const Fields kFields = Fields(Fields::Single::kData) |
                       Fields(Fields::Single::kWeights) |
                       Fields(Fields::Single::kFlags);

/**
 * Creates successive BDA buffers.
 * The resulting buffers hold values for all fields, including flags, which
 * allows testing if AppendAndWeight discards unnecessary fields.
 * @param direction_names Names of model data buffers, which should be created.
 * @param weight_value Weight value for all buffer entries.
 * @return A series of successive BDA Buffers.
 */
std::vector<std::unique_ptr<BdaBuffer>> CreateBdaBuffers(
    const std::vector<std::string>& direction_names, float weight_value) {
  std::vector<std::unique_ptr<BdaBuffer>> buffers;

  double time = kFirstTime;
  for (std::size_t i = 0; i < 3; ++i) {
    auto buffer =
        std::make_unique<BdaBuffer>(kRowsPerBuffer * kNElementsPerRow, kFields);

    // Add a non-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 0, kNChannels, kNCorrelations);
    buffer->AddRow(time + kInterval, kInterval, kInterval, 0, kNChannels,
                   kNCorrelations);

    // Add a time-averaged baseline.
    buffer->AddRow(time + kInterval / 2, kInterval * 2.0, kInterval * 2.0, 1,
                   kNChannels, kNCorrelations);

    // Add a frequency-averaged baseline.
    buffer->AddRow(time, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations);
    buffer->AddRow(time + kInterval, kInterval, kInterval, 2, kNChannels / 2,
                   kNCorrelations);

    // Fill the main data buffers with a different value per row.
    for (std::size_t row = 0; row < buffer->GetRows().size(); ++row) {
      const std::complex<float> value =
          kFirstDataValue + float(i * kRowsPerBuffer + row) * kDataIncrement;
      std::fill_n(buffer->GetData(row), buffer->GetRows()[row].GetDataSize(),
                  value);
    }

    // Set all flags to false and weights to weight_value.
    const size_t n_elements = buffer->GetNumberOfElements();
    std::fill_n(buffer->GetFlags(), n_elements, false);
    std::fill_n(buffer->GetWeights(), n_elements, weight_value);

    // Create and fill model data buffers.
    const std::complex<float>* main_data = buffer->GetData();
    for (std::size_t dir = 0; dir < direction_names.size(); ++dir) {
      const std::string& name = direction_names[dir];
      buffer->AddData(name);
      std::complex<float>* direction_data = buffer->GetData(name);
      for (std::size_t j = 0; j < n_elements; ++j) {
        direction_data[j] = main_data[j] + float(dir + 1) * kModelDataDiff;
      }
    }

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
void CheckRow(const BdaSolverBuffer& solver_buffer,
              bool has_unweighted_model_data, size_t row_index,
              const BdaBuffer::Row& reference_row,
              const std::complex<float>* reference_row_data) {
  // Get the row from the solver buffer.
  const std::vector<BdaSolverBuffer::IntervalRow>& rows =
      solver_buffer.GetIntervalRows();
  BOOST_REQUIRE_LT(row_index, rows.size());
  const BdaSolverBuffer::IntervalRow& row = rows[row_index];
  BOOST_REQUIRE(row.unweighted_data);
  if (has_unweighted_model_data) {
    BOOST_REQUIRE(!row.unweighted_model_data.empty());
    BOOST_REQUIRE(row.unweighted_model_data.front());
  } else {
    BOOST_CHECK(row.unweighted_model_data.empty());
  }
  BOOST_REQUIRE(row.weighted_data);
  BOOST_REQUIRE(!row.weighted_model_data.empty());
  BOOST_REQUIRE(row.weighted_model_data.front());

  // Compare the solver buffer row against the original row.
  CheckRowMetaData(reference_row, row);
  // Checking the value for the first correlation only should suffice.
  // Since the weights are 1 in the test, weighted data equals unweighted data.
  CheckComplex(*row.unweighted_data, *reference_row_data);
  if (has_unweighted_model_data) {
    CheckComplex(*row.unweighted_model_data.front(),
                 *reference_row_data + kModelDataDiff);
  }
  CheckComplex(*row.weighted_data, *reference_row_data);
  CheckComplex(*row.weighted_model_data.front(),
               *reference_row_data + kModelDataDiff);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bdasolverbuffer)

BOOST_AUTO_TEST_CASE(intervalrow) {
  // Test the IntervalRow constructor.
  std::complex<float> unweighted_data;
  const std::complex<float> weighted_data;
  const bool flag = false;
  const float weight = 42.0f;
  const std::size_t kNModelPointers = 42;
  std::complex<float> unweighted_model_data;
  std::complex<float> weighted_model_data;

  std::vector<std::complex<float>*> unweighted_model_ptrs(
      kNModelPointers, &unweighted_model_data);
  std::vector<const std::complex<float>*> weighted_model_ptrs(
      kNModelPointers, &weighted_model_data);
  // Check if the vectors are moved (and not copied) using their iterators,
  // which should remain equal after moving.
  const std::vector<std::complex<float>*>::iterator unweighted_model_iterator =
      unweighted_model_ptrs.begin();
  const std::vector<const std::complex<float>*>::iterator
      weighted_model_iterator = weighted_model_ptrs.begin();

  BdaSolverBuffer::IntervalRow row(
      kFirstTime, kNBaselines, kNChannels, kNCorrelations, &unweighted_data,
      &weighted_data, &flag, &weight, std::move(unweighted_model_ptrs),
      std::move(weighted_model_ptrs));
  BOOST_CHECK_EQUAL(row.time, kFirstTime);
  BOOST_CHECK_EQUAL(row.baseline_nr, kNBaselines);
  BOOST_CHECK_EQUAL(row.n_channels, kNChannels);
  BOOST_CHECK_EQUAL(row.n_correlations, kNCorrelations);
  BOOST_CHECK_EQUAL(row.unweighted_data, &unweighted_data);
  BOOST_CHECK_EQUAL(row.weighted_data, &weighted_data);
  BOOST_CHECK_EQUAL(row.flags, &flag);
  BOOST_CHECK_EQUAL(row.weights, &weight);
  BOOST_CHECK_EQUAL(row.unweighted_model_data.size(), kNModelPointers);
  BOOST_CHECK_EQUAL(row.unweighted_model_data[0], &unweighted_model_data);
  BOOST_CHECK_EQUAL(row.weighted_model_data.size(), kNModelPointers);
  BOOST_CHECK_EQUAL(row.weighted_model_data[0], &weighted_model_data);

  // Test that the original vectors were moved and not copied.
  BOOST_CHECK(unweighted_model_ptrs.begin() != unweighted_model_iterator);
  BOOST_CHECK(row.unweighted_model_data.begin() == unweighted_model_iterator);
  BOOST_CHECK(weighted_model_ptrs.begin() != weighted_model_iterator);
  BOOST_CHECK(row.weighted_model_data.begin() == weighted_model_iterator);

  // Test the move constructor of IntervalRow.
  const BdaSolverBuffer::IntervalRow moved(std::move(row));
  BOOST_CHECK_EQUAL(moved.time, kFirstTime);
  BOOST_CHECK_EQUAL(moved.baseline_nr, kNBaselines);
  BOOST_CHECK_EQUAL(moved.n_channels, kNChannels);
  BOOST_CHECK_EQUAL(moved.n_correlations, kNCorrelations);
  BOOST_CHECK_EQUAL(moved.unweighted_data, &unweighted_data);
  BOOST_CHECK_EQUAL(moved.weighted_data, &weighted_data);
  BOOST_CHECK_EQUAL(moved.flags, &flag);
  BOOST_CHECK_EQUAL(moved.weights, &weight);
  BOOST_CHECK_EQUAL(moved.unweighted_model_data.size(), kNModelPointers);
  BOOST_CHECK_EQUAL(moved.unweighted_model_data[0], &unweighted_model_data);
  BOOST_CHECK_EQUAL(moved.weighted_model_data.size(), kNModelPointers);
  BOOST_CHECK_EQUAL(moved.weighted_model_data[0], &weighted_model_data);

  // Test that the vector members were moved and not copied.
  // (When members are const, a default move constructor will copy them.)
  BOOST_CHECK(row.unweighted_model_data.begin() != unweighted_model_iterator);
  BOOST_CHECK(moved.unweighted_model_data.begin() == unweighted_model_iterator);
  BOOST_CHECK(row.weighted_model_data.begin() != weighted_model_iterator);
  BOOST_CHECK(moved.weighted_model_data.begin() == weighted_model_iterator);
}

BOOST_AUTO_TEST_CASE(constructor) {
  BdaSolverBuffer buffer(kFirstTime, kInterval, 1);
  BOOST_CHECK_EQUAL(buffer.BufferCount(), 0u);
  BOOST_CHECK(buffer.GetIntervalRows().empty());
  BOOST_CHECK(buffer.GetDone().empty());
  BOOST_CHECK_EQUAL(buffer.GetCurrentInterval(), 0);
  BOOST_CHECK_EQUAL(buffer.CurrentIntervalStart(), kFirstTime);
  BOOST_CHECK_EQUAL(buffer.IntervalDuration(), kInterval);
}

BOOST_AUTO_TEST_CASE(current_interval) {
  // Test using AdvanceInterval() on an empty buffer, and test
  // GetCurrentInterval() and CurrentIntervalStart().
  // BOOST_CHECK_CLOSE is not needed, since rounding errors should not happen.
  BdaSolverBuffer buffer(kFirstTime, kInterval, kNBaselines);
  BOOST_CHECK_EQUAL(buffer.GetCurrentInterval(), 0);
  BOOST_CHECK_EQUAL(buffer.CurrentIntervalStart(), kFirstTime);
  buffer.AdvanceInterval();
  BOOST_CHECK_EQUAL(buffer.GetCurrentInterval(), 1);
  BOOST_CHECK_EQUAL(buffer.CurrentIntervalStart(), kFirstTime + kInterval);
  buffer.AdvanceInterval();
  BOOST_CHECK_EQUAL(buffer.GetCurrentInterval(), 2);
  BOOST_CHECK_EQUAL(buffer.CurrentIntervalStart(), kFirstTime + 2 * kInterval);
}

BOOST_AUTO_TEST_CASE(append_and_weight) {
  const float kWeight = 0.8;
  const float kWeightSqrt = std::sqrt(kWeight);
  const std::vector<std::string> kDirectionNames = {"Sun42", "Moon"};
  const bool kKeepUnweightedModelData = true;

  std::vector<std::unique_ptr<BdaBuffer>> data =
      CreateBdaBuffers(kDirectionNames, kWeight);
  std::vector<BdaBuffer*> data_ptr;

  // All bda rows are in the first interval.
  BdaSolverBuffer buffer(0.0, 100.0, kNBaselines);
  for (size_t i = 0; i < data.size(); ++i) {
    data_ptr.push_back(data[i].get());
    buffer.AppendAndWeight(std::move(data[i]), kDirectionNames,
                           kKeepUnweightedModelData);
  }

  const size_t n_rows = data.size() * kRowsPerBuffer;
  BOOST_REQUIRE_EQUAL(buffer.GetIntervalRows().size(), n_rows);

  for (size_t row_index = 0; row_index < n_rows; ++row_index) {
    const BdaSolverBuffer::IntervalRow& solver_row =
        buffer.GetIntervalRows()[row_index];
    BOOST_REQUIRE(solver_row.unweighted_data);
    BOOST_REQUIRE_EQUAL(solver_row.unweighted_model_data.size(), 2);
    BOOST_REQUIRE(solver_row.unweighted_model_data[0]);
    BOOST_REQUIRE(solver_row.unweighted_model_data[1]);

    BOOST_REQUIRE(solver_row.weighted_data);
    BOOST_REQUIRE_EQUAL(solver_row.weighted_model_data.size(), 2);
    BOOST_REQUIRE(solver_row.weighted_model_data[0]);
    BOOST_REQUIRE(solver_row.weighted_model_data[1]);

    const size_t data_buffer_index = row_index / kRowsPerBuffer;
    const size_t data_row_index = row_index % kRowsPerBuffer;
    const BdaBuffer::Row& original_row =
        data_ptr[data_buffer_index]->GetRows()[data_row_index];

    CheckRowMetaData(original_row, solver_row);

    for (size_t i = 0; i < original_row.GetDataSize(); ++i) {
      const std::complex<float> inc = float(row_index) * kDataIncrement;

      const std::complex<float> unweighted_data_value = kFirstDataValue + inc;
      const std::complex<float> unweighted_model_value1 =
          unweighted_data_value + 1.0f * kModelDataDiff;
      const std::complex<float> unweighted_model_value2 =
          unweighted_data_value + 2.0f * kModelDataDiff;

      const std::complex<float> weighted_data_value =
          unweighted_data_value * kWeightSqrt;
      const std::complex<float> weighted_model_value1 =
          unweighted_model_value1 * kWeightSqrt;
      const std::complex<float> weighted_model_value2 =
          unweighted_model_value2 * kWeightSqrt;

      CheckComplex(solver_row.unweighted_data[i], unweighted_data_value);
      CheckComplex(solver_row.unweighted_model_data[0][i],
                   unweighted_model_value1);
      CheckComplex(solver_row.unweighted_model_data[1][i],
                   unweighted_model_value2);

      CheckComplex(solver_row.weighted_data[i], weighted_data_value);
      CheckComplex(solver_row.weighted_model_data[0][i], weighted_model_value1);
      CheckComplex(solver_row.weighted_model_data[1][i], weighted_model_value2);
    }
  }

  buffer.Clear();
  BOOST_TEST(buffer.GetIntervalRows().empty());
}

BOOST_AUTO_TEST_CASE(one_interval_per_buffer) {
  const float kWeight = 1.0;
  const std::vector<std::string> kDirectionNames = {"direction1"};

  std::vector<std::unique_ptr<BdaBuffer>> data =
      CreateBdaBuffers(kDirectionNames, kWeight);
  std::vector<BdaBuffer*> data_ptr;

  BdaSolverBuffer buffer(kFirstTime - kInterval / 2, kInterval * 2.0,
                         kNBaselines);
  for (size_t i = 0; i < data.size(); ++i) {
    data_ptr.push_back(data[i].get());
    buffer.AppendAndWeight(std::move(data[i]), kDirectionNames, false);
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

BOOST_DATA_TEST_CASE(buffer_with_multiple_intervals, kTrueFalseRange,
                     keep_unweighted_model_data) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;
  const std::string kModelDataName = "model";

  // Create a BDA buffer with successive rows and a single value per row.
  auto data_buffer =
      std::make_unique<BdaBuffer>(kNRows * kNCorrelations, kFields);
  data_buffer->AddData(kModelDataName);
  const BdaBuffer& data_ref = *data_buffer;

  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0, 1.0, 1.0, 1.0};
  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                        data_value.data(), nullptr, kWeight.data());
    std::copy_n(model_value.data(), model_value.size(),
                data_buffer->GetData(row, kModelDataName));
    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  BdaSolverBuffer solver_buffer(kFirstTime - kInterval / 2, kSolutionInterval,
                                1);
  solver_buffer.AppendAndWeight(std::move(data_buffer), {kModelDataName},
                                keep_unweighted_model_data);
  BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), 1u);

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, keep_unweighted_model_data, row_index,
             data_ref.GetRows()[row], data_ref.GetData(row));

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
}

BOOST_DATA_TEST_CASE(multiple_buffers_per_interval, kTrueFalseRange,
                     keep_unweighted_model_data) {
  const size_t kNRows = 42;
  const size_t kSolutionIntervalFactor = 5;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;
  const std::string kModelDataName = "model";

  // Create successive BdaBuffers with a single row.
  std::vector<BdaBuffer*> data_buffers;
  std::array<std::complex<float>, kNCorrelations> data_value{kFirstDataValue};
  std::array<std::complex<float>, kNCorrelations> model_value{kFirstDataValue +
                                                              kModelDataDiff};
  const std::array<float, kNCorrelations> kWeight{1.0, 1.0, 1.0, 1.0};

  BdaSolverBuffer solver_buffer(kFirstTime - kInterval / 2, kSolutionInterval,
                                1);

  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    auto data_buffer = std::make_unique<BdaBuffer>(kNCorrelations, kFields);
    data_buffer->AddData(kModelDataName);
    data_buffers.push_back(data_buffer.get());
    BOOST_REQUIRE(data_buffer->AddRow(time, kInterval, kInterval, 0, 1,
                                      kNCorrelations, nullptr, nullptr,
                                      kWeight.data()));
    std::copy_n(data_value.data(), kNCorrelations, data_buffer->GetData());
    std::copy_n(model_value.data(), kNCorrelations,
                data_buffer->GetData(kModelDataName));

    solver_buffer.AppendAndWeight(std::move(data_buffer), {kModelDataName},
                                  keep_unweighted_model_data);
    BOOST_CHECK_EQUAL(solver_buffer.BufferCount(), row + 1);

    time += kInterval;
    data_value[0] += kDataIncrement;
    model_value[0] += kDataIncrement;
  }

  size_t row_index = 0;
  for (size_t row = 0; row < kNRows; ++row) {
    CheckRow(solver_buffer, keep_unweighted_model_data, row_index,
             data_buffers[row]->GetRows().front(),
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

  // Test that GetDone() returns all original data buffers and that the model
  // presence corresponds to keep_unweighted_model_data.
  const std::vector<std::unique_ptr<BdaBuffer>> done = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(done.size(), data_buffers.size());
  for (size_t i = 0; i < data_buffers.size(); ++i) {
    BOOST_CHECK_EQUAL(done[i].get(), data_buffers[i]);
    BOOST_CHECK_EQUAL(done[i]->HasData(kModelDataName),
                      keep_unweighted_model_data);
  }
}

BOOST_AUTO_TEST_CASE(apply_solutions) {
  // Test full jones application on the first row,
  // test diagonal application on the second row,
  // test scalar application on the third row,
  // and test that the fourth row remains intact.
  const size_t kNRows = 4;
  const aocommon::MC2x2F kDataValue({10, 20}, {11, 21}, {12, 22}, {13, 23});
  const aocommon::MC2x2F kModelValue({1, 10}, {2, 20}, {3, 30}, {4, 40});
  // Using weights other than 1 allows verifying that solutions are applied
  // to the unweighted model and not to the weighted model.
  const std::array<float, kNCorrelations> kWeight{0.8, 1.1, 1.2, 0.9};
  const std::vector<std::vector<std::complex<double>>> kSolutions{
      {{1, 2}, {2, 2}, {3, 4}, {4, 4}, {5, 2}, {6, 2}, {7, 1}, {8, 1}}};

  const aocommon::MC2x2F kSolution1(kSolutions.front().data());
  const aocommon::MC2x2F kSolution2(kSolutions.front().data() + kNCorrelations);
  const aocommon::MC2x2F kSolvedFullJones =
      kSolution1.Multiply(kModelValue).MultiplyHerm(kSolution2);

  const aocommon::MC2x2F kDiagonalSolution1(
      std::complex<float>(kSolutions.front()[0]), 0.0f, 0.0f,
      std::complex<float>(kSolutions.front()[1]));
  const aocommon::MC2x2F kDiagonalSolution2(
      std::complex<float>(kSolutions.front()[2]), 0.0f, 0.0f,
      std::complex<float>(kSolutions.front()[3]));
  const aocommon::MC2x2F kSolvedDiagonal =
      kDiagonalSolution1.Multiply(kModelValue).MultiplyHerm(kDiagonalSolution2);

  const std::complex<float> kScalarFactor =
      std::complex<float>(kSolutions.front()[0]) *
      std::conj(std::complex<float>(kSolutions.front()[1]));
  const std::array<std::complex<float>, kNCorrelations> kSolvedScalar{
      kScalarFactor * kModelValue.Get(0), kScalarFactor * kModelValue.Get(1),
      kScalarFactor * kModelValue.Get(2), kScalarFactor * kModelValue.Get(3)};

  const std::vector<double> kStartFreqs{100, 200};
  const std::vector<std::vector<double>> kChanFreqs{{150}};
  const std::vector<int> kAnt1{0};
  const std::vector<int> kAnt2{1};

  // Create BDA input buffer with two intervals.
  const std::string kModelDataName = "model";
  auto data_buffer =
      std::make_unique<BdaBuffer>(kNRows * kNCorrelations, kFields);
  data_buffer->AddData(kModelDataName);

  double time = kFirstTime;
  for (size_t row = 0; row < kNRows; ++row) {
    data_buffer->AddRow(time, kInterval, kInterval, 0, 1, kNCorrelations,
                        nullptr, nullptr, kWeight.data());
    kDataValue.AssignTo(data_buffer->GetData(row));
    kModelValue.AssignTo(data_buffer->GetData(row, kModelDataName));
    time += kInterval;
  }

  BdaSolverBuffer solver_buffer(kFirstTime - kInterval / 2, kInterval,
                                kNBaselines);
  const bool kKeepUnweightedModelData = true;
  solver_buffer.AppendAndWeight(std::move(data_buffer), {kModelDataName},
                                kKeepUnweightedModelData);

  solver_buffer.ApplySolutions(kSolutions, kStartFreqs, 4, kAnt1, kAnt2,
                               kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.ApplySolutions(kSolutions, kStartFreqs, 2, kAnt1, kAnt2,
                               kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.ApplySolutions(kSolutions, kStartFreqs, 1, kAnt1, kAnt2,
                               kChanFreqs);
  solver_buffer.AdvanceInterval();
  solver_buffer.AdvanceInterval();

  const std::vector<std::unique_ptr<BdaBuffer>> done = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(done.size(), 1u);
  const BdaBuffer& done_buffer = *done.front();
  BOOST_REQUIRE_EQUAL(done_buffer.GetRows().size(), kNRows);
  BOOST_REQUIRE_EQUAL(done_buffer.GetDataNames().size(), 2u);
  BOOST_REQUIRE_EQUAL(done_buffer.GetDataNames()[0], "");
  BOOST_REQUIRE_EQUAL(done_buffer.GetDataNames()[1], kModelDataName);
  const std::complex<float>* model_fulljones =
      done_buffer.GetData(0, kModelDataName);
  const std::complex<float>* model_diagonal =
      done_buffer.GetData(1, kModelDataName);
  const std::complex<float>* model_scalar =
      done_buffer.GetData(2, kModelDataName);
  const std::complex<float>* model_uncorrected =
      done_buffer.GetData(3, kModelDataName);
  for (size_t c = 0; c < kNCorrelations; ++c) {
    CheckComplex(model_fulljones[c], kSolvedFullJones.Get(c));
    CheckComplex(model_diagonal[c], kSolvedDiagonal.Get(c));
    CheckComplex(model_scalar[c], kSolvedScalar[c]);
    BOOST_CHECK_EQUAL(model_uncorrected[c], kModelValue.Get(c));
  }
}

BOOST_AUTO_TEST_CASE(solution_interval_is_complete) {
  const size_t kSolutionIntervalFactor = 4;
  const double kSolutionInterval = kInterval * kSolutionIntervalFactor;
  double time = kFirstTime;

  // The first solution interval defined is [41 , 49]
  BdaSolverBuffer solver_buffer(kFirstTime - kInterval / 2, kSolutionInterval,
                                2);
  std::vector<std::complex<float>> data(kNElementsPerRow, 1.0);
  const std::vector<float> weights(kNElementsPerRow, 1.0);
  const bool flags[kNElementsPerRow]{false};

  // 1. Add intervals within the first solution interval
  auto buffer_within_first_solint =
      std::make_unique<BdaBuffer>(6 * kNElementsPerRow, kFields);
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

  solver_buffer.AppendAndWeight(std::move(buffer_within_first_solint), {},
                                false);
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 2. Add intervals on the edge of the solution interval
  auto buffer_edge_solint =
      std::make_unique<BdaBuffer>(2 * kNElementsPerRow, kFields);
  //  Baseline 0: [49, 51]
  buffer_edge_solint->AddRow(time + 4 * kInterval, kInterval, kInterval, 0,
                             kNChannels, kNCorrelations, data.data(), flags,
                             weights.data());
  //  Baseline 1: [45, 51]
  buffer_edge_solint->AddRow(time + 2.5 * kInterval, 3 * kInterval,
                             3 * kInterval, 0, kNChannels, kNCorrelations,
                             data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_edge_solint), {}, false);
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 2. Add interval completely contained in the next solution interval for
  // baseline 0
  auto buffer_within_second_solint_bl_0 =
      std::make_unique<BdaBuffer>(1 * kNElementsPerRow, kFields);
  // Baseline 0: [51, 53]
  buffer_within_second_solint_bl_0->AddRow(
      time + 5 * kInterval, kInterval, kInterval, 0, kNChannels, kNCorrelations,
      data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_within_second_solint_bl_0), {},
                                false);
  BOOST_CHECK_EQUAL(false, solver_buffer.IntervalIsComplete());

  // 3. Add interval completely contained in the next solution interval for
  // baseline 1
  auto buffer_within_second_solint_bl_1 =
      std::make_unique<BdaBuffer>(1 * kNElementsPerRow, kFields);
  // Baseline 1: [51, 53]
  buffer_within_second_solint_bl_1->AddRow(
      time + 5 * kInterval, kInterval, kInterval, 1, kNChannels, kNCorrelations,
      data.data(), flags, weights.data());

  solver_buffer.AppendAndWeight(std::move(buffer_within_second_solint_bl_1), {},
                                false);
  BOOST_CHECK_EQUAL(true, solver_buffer.IntervalIsComplete());
}

BOOST_DATA_TEST_CASE(keep_or_discard_model_data,
                     boost::unit_test::data::make({true, false}),
                     keep_model_data) {
  const std::string kSolveDirection = {"solve_direction"};
  const std::string kExtraDirection = {"extra_direction"};
  const std::complex<float> kDataValue(42.0, -13.0);
  const std::complex<float> kSolveValue = kDataValue + 1.0f * kDataIncrement;
  const std::complex<float> kExtraValue = kDataValue + 2.0f * kDataIncrement;
  const float kWeight = 0.8;  // Make AppendAndWeight apply an actual weight.

  auto buffer = std::make_unique<BdaBuffer>(kNElementsPerRow, kFields);
  const BdaBuffer& buffer_ref = *buffer;
  buffer->AddData(kSolveDirection);
  buffer->AddData(kExtraDirection);
  buffer->AddRow(kFirstTime, kInterval, kInterval, 0, kNChannels,
                 kNCorrelations);
  std::fill_n(buffer->GetData(), kNElementsPerRow, kDataValue);
  std::fill_n(buffer->GetData(kSolveDirection), kNElementsPerRow, kSolveValue);
  std::fill_n(buffer->GetData(kExtraDirection), kNElementsPerRow, kExtraValue);
  std::fill_n(buffer->GetWeights(), kNElementsPerRow, kWeight);

  BdaSolverBuffer solver_buffer(kFirstTime - kInterval / 2, kInterval, 1);
  // Only pass the solve direction to the solver buffer. It should ignore
  // the extra direction.
  BOOST_CHECK_NO_THROW(solver_buffer.AppendAndWeight(
      std::move(buffer), {kSolveDirection}, keep_model_data));

  // Before advancing, GetDone should return nothing.
  BOOST_REQUIRE(solver_buffer.GetDone().empty());

  // After advancing, GetDone() should return the original buffer.
  solver_buffer.AdvanceInterval();
  const std::vector<std::unique_ptr<BdaBuffer>> done = solver_buffer.GetDone();
  BOOST_REQUIRE_EQUAL(done.size(), 1u);
  BOOST_CHECK_EQUAL(done.front().get(), &buffer_ref);

  // Check if the buffer has the correct directions.
  const std::vector<std::string> names = buffer_ref.GetDataNames();
  if (keep_model_data) {
    // All directions should remain ('names' is sorted).
    const std::vector<std::string> kExpectedNames{"", kExtraDirection,
                                                  kSolveDirection};
    BOOST_REQUIRE_EQUAL_COLLECTIONS(names.begin(), names.end(),
                                    kExpectedNames.begin(),
                                    kExpectedNames.end());
  } else {
    // Only the solve direction should be gone.
    const std::vector<std::string> kExpectedNames{"", kExtraDirection};
    BOOST_REQUIRE_EQUAL_COLLECTIONS(names.begin(), names.end(),
                                    kExpectedNames.begin(),
                                    kExpectedNames.end());
  }

  // Check that that the data did not change.
  for (std::size_t i = 0; i < kNElementsPerRow; ++i) {
    BOOST_CHECK_EQUAL(buffer_ref.GetData()[i], kDataValue);
    BOOST_CHECK_EQUAL(buffer_ref.GetData(kExtraDirection)[i], kExtraValue);
    if (keep_model_data) {
      BOOST_CHECK_EQUAL(buffer_ref.GetData(kSolveDirection)[i], kSolveValue);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
