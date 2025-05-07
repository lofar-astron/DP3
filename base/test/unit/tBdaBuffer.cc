// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Unit tests for the BdaBuffer class.
/// @author Lars Krombeen & Maik Nijhuis

#include "tBdaBuffer.h"

#include <boost/test/unit_test.hpp>

using dp3::base::BdaBuffer;
using dp3::common::Fields;

// Implementatation of tBdaBuffer.h
namespace dp3 {
namespace base {
namespace test {

void CheckBDARowMetaData(const BdaBuffer::Row& left,
                         const BdaBuffer::Row& right) {
  BOOST_CHECK_EQUAL(left.time, right.time);
  BOOST_CHECK_EQUAL(left.interval, right.interval);
  BOOST_CHECK_EQUAL(left.exposure, right.exposure);
  BOOST_CHECK_EQUAL(left.row_nr, right.row_nr);
  BOOST_CHECK_EQUAL(left.baseline_nr, right.baseline_nr);
  BOOST_CHECK_EQUAL(left.n_channels, right.n_channels);
  BOOST_CHECK_EQUAL(left.n_correlations, right.n_correlations);
  BOOST_CHECK_EQUAL(left.GetDataSize(), right.GetDataSize());
  for (std::size_t i = 0; i < 3; ++i) {
    if (std::isnan(left.uvw[i])) {
      BOOST_CHECK(std::isnan(right.uvw[i]));
    } else {
      BOOST_CHECK_EQUAL(left.uvw[i], right.uvw[i]);
    }
  }
}

void CheckBDARowMetaData(const BdaBuffer& left, const BdaBuffer& right) {
  BOOST_REQUIRE_EQUAL(left.GetRows().size(), right.GetRows().size());
  auto left_row = left.GetRows().begin();
  auto right_row = right.GetRows().begin();
  while (left_row != left.GetRows().end()) {
    CheckBDARowMetaData(*left_row, *right_row);
    ++left_row;
    ++right_row;
  }
}

}  // namespace test
}  // namespace base
}  // namespace dp3

namespace {
const double kTime = {0.0};
const double kInterval{1.0};
const double kExposure{0.95};
const dp3::common::rownr_t kBaseRowNr{142};
const std::size_t kBaselineNr{42};
const std::size_t kNChannels{2};
const std::size_t kNCorrelations{3};
const std::size_t kDataSize{kNChannels * kNCorrelations};
const std::size_t kUnusedSpace{42};
const std::size_t k1Channel{1};
const std::size_t k1Correlation{1};
const std::size_t k1DataSize{k1Channel * k1Correlation};
const std::string kDataName = "foo";
const Fields kDataField(Fields::Single::kData);
const Fields kWeightsField(Fields::Single::kWeights);
const Fields kFlagsField(Fields::Single::kFlags);
const Fields kAllFields = kDataField | kWeightsField | kFlagsField;
const double kTimeEpsilon = 1.0e-8;
const double kHalfEpsilon = kTimeEpsilon / 2;
const double kTwoEpsilon = 2.0 * kTimeEpsilon;

void AddBasicRow(BdaBuffer& buffer, const std::string& name = "") {
  const std::complex<float> kData[kDataSize]{{1, 1}, {2, 2}, {3, 3},
                                             {4, 4}, {5, 5}, {6, 6}};
  const bool kFlags[kDataSize]{true, true, false, false, true, true};
  const float kWeights[kDataSize]{21, 22, 23, 24, 25, 26};
  const double kUvw[3]{41, 42, 43};
  buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr, kNChannels,
                kNCorrelations, kData, kFlags, kWeights, kUvw);

  if (name != "") {
    const std::size_t row = buffer.GetRows().size() - 1;
    std::copy_n(kData, kDataSize, buffer.GetData(row, name));
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(bdabuffer)

BOOST_AUTO_TEST_CASE(initialization) {
  const BdaBuffer buffer(2, kAllFields);
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 0u);
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), 2u);
  BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetRows().size(), 0u);
  BOOST_CHECK(buffer.HasData());
  const std::vector<std::string> names = buffer.GetDataNames();
  BOOST_REQUIRE_EQUAL(names.size(), 1u);
  BOOST_CHECK_EQUAL(names.front(), "");
}

BOOST_AUTO_TEST_CASE(copy_all_fields) {
  const std::complex<float> kData1[kDataSize]{{1, 1}, {2, 2}, {3, 3},
                                              {4, 4}, {5, 5}, {6, 6}};
  const std::complex<float> kData2[kDataSize]{{-1, -1}, {-2, -2}, {-3, -3},
                                              {-4, -4}, {-5, -5}, {-6, -6}};
  const bool kFlags[kDataSize]{true, false, true, true, false, true};
  const float kWeights1[kDataSize]{21, 22, 23, 24, 25, 26};
  const float kWeights2[kDataSize]{31, 32, 33, 34, 35, 36};
  const double kUvw[3]{41, 42, 43};

  BdaBuffer buffer(2 * kDataSize + kUnusedSpace, kAllFields);
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations, kData1, nullptr,
                            kWeights1, kUvw));
  buffer.SetBaseRowNr(kBaseRowNr);
  BOOST_CHECK(buffer.AddRow(kTime + 1., kInterval + 1., kExposure + 1.,
                            kBaselineNr + 1, kNChannels, kNCorrelations, kData2,
                            kFlags, kWeights2, kUvw));

  // Add an extra visibility buffer, which should also be copied.
  buffer.AddData(kDataName);
  for (std::size_t i = 0; i < 2 * kDataSize; ++i) {
    buffer.GetData(kDataName)[i] = {42.0f + i, 0.0f - i};
  }

  const BdaBuffer buffer_copy{buffer, kAllFields};

  // Verify the memory pool data in the copy.
  BOOST_CHECK_NE(buffer.GetData(), buffer_copy.GetData());
  BOOST_CHECK_NE(buffer.GetData(kDataName), buffer_copy.GetData(kDataName));
  BOOST_CHECK_NE(buffer.GetFlags(), buffer_copy.GetFlags());
  BOOST_CHECK_NE(buffer.GetWeights(), buffer_copy.GetWeights());
  BOOST_REQUIRE_NE(buffer_copy.GetData(), nullptr);
  BOOST_REQUIRE_NE(buffer_copy.GetData(kDataName), nullptr);
  BOOST_REQUIRE_NE(buffer_copy.GetFlags(), nullptr);
  BOOST_REQUIRE_NE(buffer_copy.GetWeights(), nullptr);

  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(),
                    buffer_copy.GetNumberOfElements());
  for (std::size_t i = 0; i < buffer.GetNumberOfElements(); ++i) {
    BOOST_CHECK_EQUAL(buffer.GetData()[i], buffer_copy.GetData()[i]);
    BOOST_CHECK_EQUAL(buffer.GetData(kDataName)[i],
                      buffer_copy.GetData(kDataName)[i]);
    BOOST_CHECK_EQUAL(buffer.GetFlags()[i], buffer_copy.GetFlags()[i]);
    BOOST_CHECK_EQUAL(buffer.GetWeights()[i], buffer_copy.GetWeights()[i]);
  }

  // Verify the copied rows.
  dp3::base::test::CheckBDARowMetaData(buffer, buffer_copy);

  // Verify that the original has remaining capacity, but the copy doesn't.
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), kUnusedSpace);
  BOOST_CHECK_EQUAL(buffer_copy.GetRemainingCapacity(), 0u);
}

BOOST_AUTO_TEST_CASE(copy_omit_fields) {
  BdaBuffer buffer(kDataSize + kUnusedSpace, kAllFields);
  AddBasicRow(buffer);
  buffer.SetBaseRowNr(kBaseRowNr);

  const BdaBuffer weights(buffer, kWeightsField);
  const BdaBuffer data_flags(buffer, kDataField | kFlagsField);

  // Verify the memory pool data in the copies.
  BOOST_CHECK(weights.GetData() == nullptr);
  BOOST_CHECK(weights.GetFlags() == nullptr);
  BOOST_CHECK(weights.GetWeights() != nullptr);
  BOOST_CHECK(weights.GetWeights() != buffer.GetWeights());

  BOOST_CHECK(data_flags.GetData() != nullptr);
  BOOST_CHECK(data_flags.GetFlags() != nullptr);
  BOOST_CHECK(data_flags.GetWeights() == nullptr);
  BOOST_CHECK(data_flags.GetData() != buffer.GetData());
  BOOST_CHECK(data_flags.GetFlags() != buffer.GetFlags());

  BOOST_REQUIRE_EQUAL(buffer.GetNumberOfElements(),
                      weights.GetNumberOfElements());
  BOOST_REQUIRE_EQUAL(buffer.GetNumberOfElements(),
                      data_flags.GetNumberOfElements());
  for (std::size_t i = 0; i < buffer.GetNumberOfElements(); ++i) {
    BOOST_CHECK_EQUAL(buffer.GetData()[i], data_flags.GetData()[i]);
    BOOST_CHECK_EQUAL(buffer.GetFlags()[i], data_flags.GetFlags()[i]);
    BOOST_CHECK_EQUAL(buffer.GetWeights()[i], weights.GetWeights()[i]);
  }

  // Verify the copied rows.
  dp3::base::test::CheckBDARowMetaData(buffer, weights);
  dp3::base::test::CheckBDARowMetaData(buffer, data_flags);

  // Verify that the original has remaining capacity, but the copies don't.
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), kUnusedSpace);
  BOOST_CHECK_EQUAL(weights.GetRemainingCapacity(), 0u);
  BOOST_CHECK_EQUAL(data_flags.GetRemainingCapacity(), 0u);
}

BOOST_AUTO_TEST_CASE(copy_add_weights) {
  BdaBuffer without_data_flags(kDataSize + kUnusedSpace, kWeightsField);
  AddBasicRow(without_data_flags);
  BOOST_CHECK(!without_data_flags.GetData());
  BOOST_CHECK(!without_data_flags.GetFlags());

  const BdaBuffer copy(without_data_flags, kAllFields);
  BOOST_CHECK(copy.GetData());
  BOOST_CHECK(copy.GetFlags());
  BOOST_CHECK_EQUAL(copy.GetRemainingCapacity(), 0u);
}

BOOST_AUTO_TEST_CASE(copy_add_data_flags) {
  BdaBuffer without_weights(kDataSize + kUnusedSpace, kDataField | kFlagsField);
  AddBasicRow(without_weights);
  BOOST_CHECK(!without_weights.GetWeights());

  const BdaBuffer copy(without_weights, kAllFields);
  BOOST_CHECK(copy.GetWeights());
  BOOST_CHECK_EQUAL(copy.GetRemainingCapacity(), 0u);
}

BOOST_AUTO_TEST_CASE(add_row_all_fields) {
  const std::complex<float> kData1[kDataSize]{{1, 1}, {2, 2}, {3, 3},
                                              {4, 4}, {5, 5}, {6, 6}};
  const std::complex<float> kData2[k1DataSize]{{7, 7}};
  const bool kFlags1[kDataSize]{true, false, true, true, false, true};
  const bool kFlags2[k1DataSize]{false};
  const float kWeights1[kDataSize]{21, 22, 23, 24, 25, 26};
  const float kWeights2[k1DataSize]{27};
  const double kUvw1[3]{31, 32, 33};
  const double kUvw2[3]{41, 42, 43};

  BdaBuffer buffer(kDataSize + k1DataSize, kAllFields);

  // Add rows and verify that the capacity is reached at the third row.
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations, kData1, kFlags1,
                            kWeights1, kUvw1));
  buffer.SetBaseRowNr(kBaseRowNr);
  BOOST_CHECK(buffer.AddRow(kTime + 1.0, kInterval + 1.0, kExposure + 1.0,
                            kBaselineNr + 1, k1Channel, k1Correlation, kData2,
                            kFlags2, kWeights2, kUvw2));
  BOOST_CHECK(!buffer.AddRow(kTime + 2.0, kInterval, kExposure, kBaselineNr + 2,
                             k1Channel, k1Correlation));

  // Verify the data in the buffer.
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), kDataSize + k1DataSize);
  BOOST_CHECK(buffer.GetData());
  BOOST_CHECK(buffer.GetFlags());
  BOOST_CHECK(buffer.GetWeights());

  BOOST_CHECK_EQUAL(buffer.GetData(0), buffer.GetData());
  BOOST_CHECK_EQUAL(buffer.GetData(1), buffer.GetData() + kDataSize);
  BOOST_CHECK_EQUAL(buffer.GetFlags(0), buffer.GetFlags());
  BOOST_CHECK_EQUAL(buffer.GetFlags(1), buffer.GetFlags() + kDataSize);
  BOOST_CHECK_EQUAL(buffer.GetWeights(0), buffer.GetWeights());
  BOOST_CHECK_EQUAL(buffer.GetWeights(1), buffer.GetWeights() + kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i) {
    BOOST_CHECK_EQUAL(buffer.GetData()[i], kData1[i]);
    BOOST_CHECK_EQUAL(buffer.GetFlags()[i], kFlags1[i]);
    BOOST_CHECK_EQUAL(buffer.GetWeights()[i], kWeights1[i]);
  }
  for (std::size_t i = 0; i < k1DataSize; ++i) {
    BOOST_CHECK_EQUAL(buffer.GetData(1)[i], kData2[i]);
    BOOST_CHECK_EQUAL(buffer.GetFlags(1)[i], kFlags2[i]);
    BOOST_CHECK_EQUAL(buffer.GetWeights(1)[i], kWeights2[i]);
  }

  // Verify the rows.
  const auto& rows = buffer.GetRows();
  BOOST_CHECK_EQUAL(rows.size(), 2u);
  BOOST_CHECK_EQUAL(rows[0].time, kTime);
  BOOST_CHECK_EQUAL(rows[0].interval, kInterval);
  BOOST_CHECK_EQUAL(rows[0].row_nr, kBaseRowNr);
  BOOST_CHECK_EQUAL(rows[0].baseline_nr, kBaselineNr);
  BOOST_CHECK_EQUAL(rows[0].n_channels, kNChannels);
  BOOST_CHECK_EQUAL(rows[0].n_correlations, kNCorrelations);
  BOOST_CHECK_EQUAL(rows[0].GetDataSize(), kDataSize);

  BOOST_CHECK_EQUAL(rows[1].time, kTime + 1.);
  BOOST_CHECK_EQUAL(rows[1].interval, kInterval + 1.);
  BOOST_CHECK_EQUAL(rows[1].row_nr, kBaseRowNr + 1);
  BOOST_CHECK_EQUAL(rows[1].baseline_nr, kBaselineNr + 1);
  BOOST_CHECK_EQUAL(rows[1].n_channels, k1Channel);
  BOOST_CHECK_EQUAL(rows[1].n_correlations, k1Correlation);
  BOOST_CHECK_EQUAL(rows[1].GetDataSize(), k1DataSize);

  BOOST_CHECK_EQUAL(rows[0].offset, 0);
  BOOST_CHECK_EQUAL(rows[1].offset, buffer.GetData(1) - buffer.GetData());
  BOOST_CHECK_EQUAL(rows[1].offset, buffer.GetFlags(1) - buffer.GetFlags());
  BOOST_CHECK_EQUAL(rows[1].offset, buffer.GetWeights(1) - buffer.GetWeights());

  for (std::size_t i = 0; i < 3; ++i) {
    BOOST_CHECK_EQUAL(rows[0].uvw[i], kUvw1[i]);
    BOOST_CHECK_EQUAL(rows[1].uvw[i], kUvw2[i]);
  }
}

BOOST_AUTO_TEST_CASE(add_row_no_fields) {
  BdaBuffer buffer(kDataSize, kAllFields);
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations));
  buffer.SetBaseRowNr(kBaseRowNr);
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), 0u);
  BOOST_CHECK(buffer.GetData());
  BOOST_CHECK(buffer.GetFlags());
  BOOST_CHECK(buffer.GetWeights());

  // Check that the memory pools hold default values.
  for (std::size_t i = 0; i < kDataSize; ++i) {
    BOOST_CHECK_EQUAL(buffer.GetData()[i], std::complex<float>(0.0, 0.0));
    BOOST_CHECK_EQUAL(buffer.GetFlags()[i], false);
    BOOST_CHECK_EQUAL(buffer.GetWeights()[i], 0.0);
  }

  // Verify the row.
  const auto& rows = buffer.GetRows();
  BOOST_CHECK_EQUAL(rows.size(), 1u);
  const auto& row = rows.front();
  BOOST_CHECK_EQUAL(row.GetDataSize(), kDataSize);
  BOOST_CHECK_EQUAL(row.time, kTime);
  BOOST_CHECK_EQUAL(row.interval, kInterval);
  BOOST_CHECK_EQUAL(row.row_nr, kBaseRowNr);
  BOOST_CHECK_EQUAL(row.baseline_nr, kBaselineNr);
  BOOST_CHECK_EQUAL(row.n_channels, kNChannels);
  BOOST_CHECK_EQUAL(row.n_correlations, kNCorrelations);
  BOOST_CHECK_EQUAL(row.offset, 0);
  for (std::size_t i = 0; i < 3; ++i) {
    BOOST_CHECK(std::isnan(row.uvw[i]));
  }
}

BOOST_AUTO_TEST_CASE(disabled_fields) {
  BdaBuffer buffer(kDataSize, Fields());
  AddBasicRow(buffer);

  // Verify the data pointers from buffer.Get*().
  BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);

  BOOST_CHECK_EQUAL(buffer.GetData(0), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetFlags(0), nullptr);
  BOOST_CHECK_EQUAL(buffer.GetWeights(0), nullptr);

  // Verify functions related to multiple visibility buffers.
  BOOST_CHECK(!buffer.HasData());
  BOOST_CHECK(buffer.GetDataNames().empty());

  // Verify the offset in the row.
  const auto& rows = buffer.GetRows();
  BOOST_REQUIRE_EQUAL(rows.size(), 1u);
  const auto& row = rows.front();
  BOOST_CHECK_EQUAL(row.offset, 0);
}

BOOST_AUTO_TEST_CASE(add_wrong_ordering) {
  const double kTime1 = 4.5;
  const double kInterval1 = 3.0;
  const double kTime2 = 1.5;
  const double kInterval2 = 3.0;

  BdaBuffer buffer(2 * kDataSize, kAllFields);

  BOOST_CHECK(buffer.AddRow(kTime1, kInterval1, kInterval1, kBaselineNr,
                            kNChannels, kNCorrelations));
  // end time equals start time of last row -> invalid ordering.
  BOOST_CHECK_THROW(buffer.AddRow(kTime2, kInterval2, kInterval2, kBaselineNr,
                                  kNChannels, kNCorrelations),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(add_overlap) {
  const double kTime1 = 3.0;
  const double kInterval1 = 3.0;
  const double kTimePartialOverlap = 3.5;
  const double kIntervalPartialOverlap = 2.0;
  const double kTimeFullOverlap = 1.0;
  const double kIntervalFullOverlap = 10.0;

  BdaBuffer buffer(3 * kDataSize, kAllFields);

  BOOST_CHECK(buffer.AddRow(kTime1, kInterval1, kInterval1, kBaselineNr,
                            kNChannels, kNCorrelations));
  // Add partially overlapping row, starting before the last row.
  BOOST_CHECK(buffer.AddRow(kTimePartialOverlap, kIntervalPartialOverlap,
                            kIntervalPartialOverlap, kBaselineNr, kNChannels,
                            kNCorrelations));
  // Add fully overlapping row, starting before the last row.
  BOOST_CHECK(buffer.AddRow(kTimeFullOverlap, kIntervalFullOverlap,
                            kIntervalFullOverlap, kBaselineNr, kNChannels,
                            kNCorrelations));
}

BOOST_AUTO_TEST_CASE(add_main_data) {
  // Reference scenerio: Add a row to a buffer without visibilities.
  BdaBuffer without_data(kDataSize, Fields());
  AddBasicRow(without_data);
  BOOST_CHECK_EQUAL(without_data.GetRows().size(), 1u);
  BOOST_CHECK_EQUAL(without_data.GetData(), nullptr);

  // Test scenerio: Call AddData() before adding the row.
  BdaBuffer with_data(kDataSize, Fields());
  with_data.AddData();
  AddBasicRow(with_data);
  BOOST_CHECK_EQUAL(with_data.GetRows().size(), 1u);
  BOOST_REQUIRE_NE(with_data.GetData(), nullptr);
  for (std::size_t i = 0; i < kDataSize; ++i) {
    BOOST_CHECK_EQUAL(with_data.GetData()[i],
                      std::complex<float>(i + 1, i + 1));
  }
}

BOOST_AUTO_TEST_CASE(add_data_before_adding_rows) {
  BdaBuffer buffer(2 * kDataSize, kAllFields);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kDataName));
  BOOST_CHECK(!buffer.GetData(kDataName));

  buffer.AddData(kDataName);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(buffer.HasData(kDataName));
  // With no rows added, all data buffers are still empty.
  BOOST_CHECK(!buffer.GetData());
  BOOST_CHECK(!buffer.GetData(kDataName));

  AddBasicRow(buffer, kDataName);
  BOOST_CHECK(buffer.GetData());
  BOOST_CHECK(buffer.GetData(kDataName));
  BOOST_CHECK_NE(buffer.GetData(kDataName), buffer.GetData());
  BOOST_CHECK_EQUAL(buffer.GetData(0, kDataName), buffer.GetData(kDataName));

  AddBasicRow(buffer, kDataName);
  BOOST_CHECK_EQUAL(buffer.GetData(1, kDataName),
                    buffer.GetData(kDataName) + kDataSize);

  // Check the values in the main and new data buffers.
  for (std::size_t row = 0; row < 2; ++row) {
    for (std::size_t i = 0; i < kDataSize; ++i) {
      const std::complex<float> expected_value(i + 1, i + 1);
      BOOST_CHECK_EQUAL(buffer.GetData(row)[i], expected_value);
      BOOST_CHECK_EQUAL(buffer.GetData(row, kDataName)[i], expected_value);
    }
  }
}

BOOST_AUTO_TEST_CASE(add_data_after_adding_rows) {
  BdaBuffer buffer(2 * kDataSize, kAllFields);
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 1,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kDataName));
  BOOST_CHECK(!buffer.GetData(kDataName));

  buffer.AddData(kDataName);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(buffer.HasData(kDataName));
  BOOST_CHECK(buffer.GetData(kDataName));
  BOOST_CHECK_EQUAL(buffer.GetData(0, kDataName), buffer.GetData(kDataName));
  BOOST_CHECK_EQUAL(buffer.GetData(1, kDataName),
                    buffer.GetData(kDataName) + kDataSize);

  // Checking visibility values is not possible, since AddData does not
  // initialize values.
}

/// Fixture for move tests, which defines source and target buffers with data.
struct SourceTargetFixture {
  SourceTargetFixture()
      : source(kDataSize, kDataField), target(kDataSize, kDataField) {
    AddBasicRow(source);
    AddBasicRow(target);
  }

  BdaBuffer source;
  BdaBuffer target;
};

BOOST_FIXTURE_TEST_CASE(move_main_data, SourceTargetFixture) {
  const BdaBuffer source_copy(source, kDataField);

  std::fill_n(target.GetData(), kDataSize, std::complex<float>(0.0, 0.0));

  target.MoveData(source);

  BOOST_CHECK_EQUAL_COLLECTIONS(source_copy.GetData(),
                                source_copy.GetData() + kDataSize,
                                target.GetData(), target.GetData() + kDataSize);
  BOOST_CHECK(source.GetDataNames().empty());
}

BOOST_FIXTURE_TEST_CASE(move_to_extra_data, SourceTargetFixture) {
  const std::string kDataName = "test_model_data";

  target.MoveData(source, "", kDataName);

  BOOST_CHECK_EQUAL_COLLECTIONS(target.GetData(), target.GetData() + kDataSize,
                                target.GetData(kDataName),
                                target.GetData(kDataName) + kDataSize);
  BOOST_CHECK(source.GetDataNames().empty());
}

BOOST_FIXTURE_TEST_CASE(move_within_buffer, SourceTargetFixture) {
  const std::string kDataName = "test_model_data";

  source.MoveData(source);  // Move main buffer to itself.
  BOOST_CHECK(source.HasData());

  source.MoveData(source, "", kDataName);  // Move main buffer to extra buffer.
  BOOST_CHECK(!source.HasData());
  BOOST_CHECK(source.HasData(kDataName));

  source.MoveData(source, kDataName, kDataName);  // Move extra buf to itself.
  BOOST_CHECK(!source.HasData());
  BOOST_CHECK(source.HasData(kDataName));

  source.MoveData(source, kDataName, "");  // Move extra buffer to main buffer.
  BOOST_CHECK(source.HasData());
  BOOST_CHECK(!source.HasData(kDataName));

  // Check if the data is still the same after all moves.
  // This test only changes 'source' so it can compare against 'target'.
  BOOST_CHECK_EQUAL_COLLECTIONS(target.GetData(), target.GetData() + kDataSize,
                                source.GetData(), source.GetData() + kDataSize);
}

BOOST_AUTO_TEST_CASE(remove_main_data) {
  BdaBuffer buffer(kDataSize, kDataField);
  AddBasicRow(buffer);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(buffer.GetData());
  BOOST_CHECK(buffer.GetDataNames().size() == 1u &&
              buffer.GetDataNames().front() == "");

  buffer.RemoveData();
  BOOST_CHECK(!buffer.HasData());
  BOOST_CHECK(!buffer.GetData());
  BOOST_CHECK(buffer.GetDataNames().empty());

  BOOST_CHECK_NO_THROW(buffer.RemoveData());
}

BOOST_AUTO_TEST_CASE(remove_named_data) {
  BdaBuffer buffer(kDataSize, Fields());
  buffer.AddData(kDataName);
  AddBasicRow(buffer);
  BOOST_CHECK(!buffer.HasData());
  BOOST_CHECK(buffer.HasData(kDataName));
  BOOST_CHECK(buffer.GetDataNames().size() == 1u &&
              buffer.GetDataNames().front() == kDataName);
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), kDataSize);

  buffer.RemoveData(kDataName);
  BOOST_CHECK(!buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kDataName));
  BOOST_CHECK(buffer.GetDataNames().empty());
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), kDataSize);

  BOOST_CHECK_NO_THROW(buffer.RemoveData(kDataName));
}

BOOST_AUTO_TEST_CASE(clear) {
  const auto check_data_names = [](const BdaBuffer& buffer) {
    BOOST_CHECK(buffer.HasData());
    BOOST_CHECK(buffer.HasData(kDataName));

    const std::vector<std::string> expected_names = {"", kDataName};
    const std::vector<std::string> names = buffer.GetDataNames();
    BOOST_CHECK_EQUAL_COLLECTIONS(names.begin(), names.end(),
                                  expected_names.begin(), expected_names.end());
  };

  BdaBuffer buffer(3 * kDataSize, kAllFields);
  buffer.AddData(kDataName);

  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 1,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 2,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(!buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 3,
                             kNChannels, kNCorrelations));
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 3 * kDataSize);
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), 0u);
  BOOST_CHECK_EQUAL(buffer.GetRows().size(), 3u);
  check_data_names(buffer);

  buffer.Clear();
  BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 0u);
  BOOST_CHECK_EQUAL(buffer.GetRemainingCapacity(), 3 * kDataSize);
  BOOST_CHECK_EQUAL(buffer.GetRows().size(), 0u);
  // Clearing should only clear the content, not the data buffers themselves.
  check_data_names(buffer);

  // Check that 3 rows can be added again.
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 1,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 2,
                            kNChannels, kNCorrelations));
  BOOST_CHECK(!buffer.AddRow(kTime, kInterval, kExposure, kBaselineNr + 3,
                             kNChannels, kNCorrelations));
}

BOOST_AUTO_TEST_CASE(time_is_less) {
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLess(5., 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLess(15., 10.), false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLess(10. - kTwoEpsilon, 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLess(10. - kHalfEpsilon, 10.), false);
}

BOOST_AUTO_TEST_CASE(time_is_less_equal) {
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLessEqual(5., 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLessEqual(15., 10.), false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLessEqual(10. + kTwoEpsilon, 10.), false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsLessEqual(10. + kHalfEpsilon, 10.), true);
}

BOOST_AUTO_TEST_CASE(time_is_greater_equal) {
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsGreaterEqual(5., 10.), false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsGreaterEqual(15., 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsGreaterEqual(10., 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsGreaterEqual(10. - kTwoEpsilon, 10.),
                    false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsGreaterEqual(10. - kHalfEpsilon, 10.),
                    true);
}

BOOST_AUTO_TEST_CASE(time_is_equal) {
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsEqual(10., 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsEqual(20., 10.), false);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsEqual(10. + kHalfEpsilon, 10.), true);
  BOOST_CHECK_EQUAL(BdaBuffer::TimeIsEqual(10. - kTwoEpsilon, 10.), false);
}

BOOST_AUTO_TEST_CASE(metadata_comparison) {
  const std::size_t kPoolSize = kBaselineNr * kNChannels * kNCorrelations;

  BdaBuffer buffer(kPoolSize, Fields());
  buffer.AddRow(1.0, 2.0, 2.0, 1, kNChannels, kNCorrelations);
  buffer.AddRow(1.0, 2.0, 2.0, 2, kNChannels, kNCorrelations);

  BdaBuffer buffer_different_n_rows(kPoolSize, Fields());
  buffer_different_n_rows.AddRow(1.0, 2.0, 2.0, 1, kNChannels, kNCorrelations);

  BdaBuffer buffer_different_baseline_ordering(
      kBaselineNr * kNChannels * kNCorrelations, Fields());
  buffer_different_baseline_ordering.AddRow(1.0, 2.0, 2.0, 2, kNChannels,
                                            kNCorrelations);
  buffer_different_baseline_ordering.AddRow(1.0, 2.0, 2.0, 1, kNChannels,
                                            kNCorrelations);

  BdaBuffer buffer_different_time(kPoolSize, Fields());
  buffer_different_time.AddRow(2.0, 2.0, 2.0, 2, kNChannels, kNCorrelations);
  buffer_different_time.AddRow(2.0, 2.0, 2.0, 1, kNChannels, kNCorrelations);

  BdaBuffer buffer_different_interval(kPoolSize, Fields());
  buffer_different_interval.AddRow(1.0, 3.0, 2.0, 2, kNChannels,
                                   kNCorrelations);
  buffer_different_interval.AddRow(1.0, 3.0, 2.0, 1, kNChannels,
                                   kNCorrelations);

  BdaBuffer buffer_different_exposure(kPoolSize, Fields());
  buffer_different_exposure.AddRow(1.0, 2.0, 3.0, 2, kNChannels,
                                   kNCorrelations);
  buffer_different_exposure.AddRow(1.0, 2.0, 3.0, 1, kNChannels,
                                   kNCorrelations);

  BdaBuffer buffer_different_nchan(
      kBaselineNr * (kNChannels + 1) * kNCorrelations, Fields());
  buffer_different_nchan.AddRow(1.0, 2.0, 2.0, 2, kNChannels + 1,
                                kNCorrelations);
  buffer_different_nchan.AddRow(1.0, 2.0, 2.0, 1, kNChannels + 1,
                                kNCorrelations);

  BdaBuffer buffer_different_ncorr(
      kBaselineNr * kNChannels * (kNCorrelations + 1), Fields());
  buffer_different_ncorr.AddRow(1.0, 2.0, 2.0, 2, kNChannels,
                                kNCorrelations + 1);
  buffer_different_ncorr.AddRow(1.0, 2.0, 2.0, 1, kNChannels,
                                kNCorrelations + 1);

  BdaBuffer buffer_equal(kPoolSize, Fields());
  buffer_equal.AddRow(1.0, 2.0, 2.0, 1, kNChannels, kNCorrelations);
  buffer_equal.AddRow(1.0, 2.0, 2.0, 2, kNChannels, kNCorrelations);

  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_n_rows));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_baseline_ordering));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_time));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_interval));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_exposure));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_nchan));
  BOOST_CHECK(!buffer.IsMetadataEqual(buffer_different_ncorr));
  BOOST_CHECK(buffer.IsMetadataEqual(buffer_equal));
}

BOOST_AUTO_TEST_SUITE_END()
