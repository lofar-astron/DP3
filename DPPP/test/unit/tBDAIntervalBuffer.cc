// Copyexpected (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Unit tests for the BDAIntervalBuffer class.
/// @author Maik Nijhuis

#include "../../BDAIntervalBuffer.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

using DP3::DPPP::BDABuffer;
using DP3::DPPP::BDAIntervalBuffer;

namespace {
  const double kTime = 0.0;
  const double kInterval = 10.0;
  const double kRowInterval = 2.0;
  const double kMaxRowInterval = 20.0;
  const float kWeight = 2.0;

  const DP3::rownr_t kRowNr = 42;
  const std::size_t kBaselineNr = 42;
  const std::size_t kNChannels = 3;
  const std::size_t kNCorrelations = 2;
  const std::size_t kDataSize = kNChannels * kNCorrelations;

  /**
   * Add a row and fill its data values.
   */
  void AddRow(BDABuffer& buffer, double time, double interval,
              DP3::rownr_t row_nr, std::size_t baseline_nr,
              bool flag, float weight)
  {
    const bool flags[kDataSize]{ flag };

    BOOST_CHECK(buffer.AddRow(time, interval, row_nr, baseline_nr,
                              kNChannels, kNCorrelations,
                              nullptr, flags, nullptr, flags));

    const BDABuffer::Row& row = buffer.GetRows().back();
    for (std::size_t i = 0; i < row.GetDataSize(); ++i) {
      const float value = row_nr * kDataSize + i;
      row.data_[i] = { value, -value };
      row.weights_[i] = weight;
    }
  }

  /**
   * Add a 'short' row, with interval 2, to a BDABuffer.
   * @param weight_factor The contribution factor for the row.
   */
  void AddShortRow(BDABuffer& buffer, double time, DP3::rownr_t row_nr,
                   bool flag, float weight_factor = 1.0f)
  {
    AddRow(buffer, time, kRowInterval * weight_factor, row_nr, kBaselineNr,
           flag, kWeight * weight_factor);
  }

  /**
   * Add a 'long' row, with interval 4, to a BDABuffer.
   * @param weight_factor The contribution factor for the row.
   */
  void AddLongRow(BDABuffer& buffer, double time, DP3::rownr_t row_nr,
                  float weight_factor = 1.0f)
  {
    AddRow(buffer, time, 2.0 * kRowInterval * weight_factor, row_nr,
           kBaselineNr + 1, false, kWeight * weight_factor);
  }

  /**
   * Advance the interval buffer and check that GetBuffer() returns
   * an empty buffer.
   */
  void CheckAdvanceBeyondEnd(BDAIntervalBuffer& interval)
  {
    // Advance the interval: GetBuffer() should return an empty buffer now.
    interval.Advance(kInterval);
    BOOST_CHECK(!interval.IsComplete());
    BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 0u);
    std::unique_ptr<BDABuffer> buffer = interval.GetBuffer();
    BOOST_CHECK(buffer && buffer->GetRows().empty());
  }

  void CheckBuffer(const BDAIntervalBuffer& interval, const BDABuffer& expected)
  {
    const std::unique_ptr<BDABuffer> result = interval.GetBuffer();

    BOOST_CHECK_EQUAL(result->GetRows().size(), expected.GetRows().size());
    BOOST_CHECK_EQUAL(result->GetNumberOfElements(), expected.GetNumberOfElements());

    auto result_row_it = result->GetRows().begin();
    auto expected_row_it = expected.GetRows().begin();
    while (result_row_it != result->GetRows().end()) {
      BOOST_CHECK_EQUAL(result_row_it->time_, expected_row_it->time_);
      BOOST_CHECK_EQUAL(result_row_it->interval_, expected_row_it->interval_);
      BOOST_CHECK_EQUAL(result_row_it->row_nr_, expected_row_it->row_nr_);
      BOOST_CHECK_EQUAL(result_row_it->baseline_nr_, expected_row_it->baseline_nr_);
      // Ignore other row members in this test: The checked values already clearly
      // indicate that the intervalbuffer copied the correct row into the result.

      for (std::size_t i = 0; i < kDataSize; ++i) {
        // Data should be absolutely equal: The intervalbuffer should not
        // perform any computations on it and only copy the data.
        BOOST_CHECK_EQUAL(result_row_it->data_[i], expected_row_it->data_[i]);
        BOOST_CHECK_EQUAL(result_row_it->flags_[i], expected_row_it->flags_[i]);
        // Weights may differ slightly because of rounding errors.
        BOOST_CHECK(BDABuffer::TimeIsEqual(result_row_it->weights_[i],
                                           expected_row_it->weights_[i]));
        BOOST_CHECK_EQUAL(result_row_it->full_res_flags_[i],
                          expected_row_it->full_res_flags_[i]);
      }

      ++result_row_it;
      ++expected_row_it;
    }

    // Verify that the result has no remaining capacacity.
    BOOST_CHECK(!result->AddRow(kTime + 42.0, kRowInterval, 42, kBaselineNr, 1, 1));
  }

  /**
   * Common part of the add_cross_boundary tests, for checking the result.
   * @param interval The interval that is under test.
   */
  void CheckCrossBoundary(BDAIntervalBuffer& interval)
  {
    auto expected_buffer1 = boost::make_unique<BDABuffer>(kDataSize * 8);
    AddShortRow(*expected_buffer1, kTime + 1.0, kRowNr + 0, false);
    AddShortRow(*expected_buffer1, kTime + 3.0, kRowNr + 1, true);
    AddLongRow(*expected_buffer1, kTime + 1.0, kRowNr + 2);
    AddShortRow(*expected_buffer1, kTime + 5.0, kRowNr + 3, false);
    AddShortRow(*expected_buffer1, kTime + 7.0, kRowNr + 4, true);
    AddLongRow(*expected_buffer1, kTime + 5.0, kRowNr + 5);
    // This short row has 1/2 in the first interval.
    AddShortRow(*expected_buffer1, kTime + 9.0, kRowNr + 6, false, 1.0f / 2.0f);
    // This long row has 1/4 in the first interval.
    AddLongRow(*expected_buffer1, kTime + 9.0, kRowNr + 8, 1.0f / 4.0f);

    auto expected_buffer2 = boost::make_unique<BDABuffer>(kDataSize * 6);
    // This short row has 1/2 in the second interval.
    AddShortRow(*expected_buffer2, kTime + 10.0, kRowNr + 6, false, 1.0f / 2.0f);
    AddShortRow(*expected_buffer2, kTime + 11.0, kRowNr + 7, true);
    // This long row has 3/4 in the second interval.
    AddLongRow(*expected_buffer2, kTime + 10.0, kRowNr + 8, 3.0f / 4.0f);
    AddShortRow(*expected_buffer2, kTime + 13.0, kRowNr + 9, false);
    AddShortRow(*expected_buffer2, kTime + 15.0, kRowNr + 10, true);
    AddLongRow(*expected_buffer2, kTime + 13.0, kRowNr + 11);

    CheckBuffer(interval, *expected_buffer1);

    interval.Advance(kInterval);

    CheckBuffer(interval, *expected_buffer2);

    CheckAdvanceBeyondEnd(interval);
  }
}

BOOST_AUTO_TEST_SUITE( bdaintervalbuffer )

BOOST_AUTO_TEST_CASE( initialization )
{
  const BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 0u);
  std::unique_ptr<BDABuffer> resultBuffer = interval.GetBuffer();
  BOOST_CHECK(resultBuffer && resultBuffer->GetRows().empty());
}

BOOST_AUTO_TEST_CASE( add_single )
{
  BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};

  auto buffer = boost::make_unique<BDABuffer>(kDataSize);
  AddShortRow(*buffer, kTime, 42, false);

  interval.AddBuffer(*buffer);
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 1u);
  CheckBuffer(interval, *buffer);
}

BOOST_AUTO_TEST_CASE( add_multiple )
{
  BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};

  auto buffer1 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer2 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer3 = boost::make_unique<BDABuffer>(kDataSize);
  auto combined = boost::make_unique<BDABuffer>(kDataSize * 3);

  AddShortRow(*buffer1, kTime + 1.0, 41, false);
  AddShortRow(*buffer2, kTime + 2.0, 42, false);
  AddShortRow(*buffer3, kTime + 3.0, 43, false);

  AddShortRow(*combined, kTime + 1.0, 41, false);
  AddShortRow(*combined, kTime + 2.0, 42, false);
  AddShortRow(*combined, kTime + 3.0, 43, false);

  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer1));
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer2));
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer3));
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 3u);
  CheckBuffer(interval, *combined);
}

BOOST_AUTO_TEST_CASE( add_cross_boundary_single )
{
  BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};

  auto buffer = boost::make_unique<BDABuffer>(kDataSize * 12);
  AddShortRow(*buffer, kTime + 1.0, kRowNr + 0, false);
  AddShortRow(*buffer, kTime + 3.0, kRowNr + 1, true);
  AddLongRow(*buffer, kTime + 1.0, kRowNr + 2);
  AddShortRow(*buffer, kTime + 5.0, kRowNr + 3, false);
  AddShortRow(*buffer, kTime + 7.0, kRowNr + 4, true);
  AddLongRow(*buffer, kTime + 5.0, kRowNr + 5);
  AddShortRow(*buffer, kTime + 9.0, kRowNr + 6, false);
  AddShortRow(*buffer, kTime + 11.0, kRowNr + 7, true);
  AddLongRow(*buffer, kTime + 9.0, kRowNr + 8);
  AddShortRow(*buffer, kTime + 13.0, kRowNr + 9, false);
  AddShortRow(*buffer, kTime + 15.0, kRowNr + 10, true);
  AddLongRow(*buffer, kTime + 13.0, kRowNr + 11);

  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer));
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 1u);
  CheckCrossBoundary(interval);
}

BOOST_AUTO_TEST_CASE( add_cross_boundary_multiple )
{
  BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};

  auto buffer1 = boost::make_unique<BDABuffer>(kDataSize * 3);
  AddShortRow(*buffer1, kTime + 1.0, kRowNr + 0, false);
  AddShortRow(*buffer1, kTime + 3.0, kRowNr + 1, true);
  AddLongRow(*buffer1, kTime + 1.0, kRowNr + 2);
  auto buffer2 = boost::make_unique<BDABuffer>(kDataSize * 3);
  AddShortRow(*buffer2, kTime + 5.0, kRowNr + 3, false);
  AddShortRow(*buffer2, kTime + 7.0, kRowNr + 4, true);
  AddLongRow(*buffer2, kTime + 5.0, kRowNr + 5);
  auto buffer3 = boost::make_unique<BDABuffer>(kDataSize * 3);
  AddShortRow(*buffer3, kTime + 9.0, kRowNr + 6, false);
  AddShortRow(*buffer3, kTime + 11.0, kRowNr + 7, true);
  AddLongRow(*buffer3, kTime + 9.0, kRowNr + 8);
  auto buffer4 = boost::make_unique<BDABuffer>(kDataSize * 3);
  AddShortRow(*buffer4, kTime + 13.0, kRowNr + 9, false);
  AddShortRow(*buffer4, kTime + 15.0, kRowNr + 10, true);
  AddLongRow(*buffer4, kTime + 13.0, kRowNr + 11);

  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer1));
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer2));
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer3));
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer4));
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 4u);
  CheckCrossBoundary(interval);
}

BOOST_AUTO_TEST_CASE( add_longer_than_interval )
{
  const double kRowInterval = 16.0;
  const double kRowInterval1 = 1.0; // Overlap in the first BDA interval.
  const double kRowInterval2 = kInterval; // Overlap in the second BDA interval.
  const double kRowInterval3 = 5.0; // Overlap in the third BDA interval.

  BDAIntervalBuffer interval{0.0, kInterval, kMaxRowInterval};

  auto buffer = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*buffer, 9.0, kRowInterval, kRowNr, kBaselineNr, true, kWeight);

  auto expected1 = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*expected1, 9.0, kRowInterval1, kRowNr, kBaselineNr,
         true, kWeight * kRowInterval1 / kRowInterval);

  // Rows that fully overlap the interval do not have an adjusted weight.
  auto expected2 = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*expected2, kInterval, kRowInterval2, kRowNr, kBaselineNr, true, kWeight);

  auto expected3 = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*expected3, 2.0 * kInterval, kRowInterval3, kRowNr,
         kBaselineNr, true, kWeight * kRowInterval3 / kRowInterval);


  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer));
  CheckBuffer(interval, *expected1);

  interval.Advance(kInterval);
  CheckBuffer(interval, *expected2);

  interval.Advance(kInterval);
  CheckBuffer(interval, *expected3);
}

BOOST_AUTO_TEST_CASE( add_max_interval )
{
  auto buffer_max = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*buffer_max, kTime, kMaxRowInterval, kRowNr, kBaselineNr, true, kWeight);

  auto buffer_over_max = boost::make_unique<BDABuffer>(kDataSize);
  AddRow(*buffer_over_max,
         kTime + kMaxRowInterval + 4.0, // Do not overlap the row in buffer_max.
         kMaxRowInterval + 0.01, // Set a too high interval.
         kRowNr, kBaselineNr, true, kWeight);


  BDAIntervalBuffer interval{0.0, kInterval, kMaxRowInterval};
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer_max));
  BOOST_CHECK_THROW(interval.AddBuffer(*buffer_over_max), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( is_complete )
{
  BDAIntervalBuffer interval{kTime, kInterval, kMaxRowInterval};

  auto buffer1 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer2 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer3 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer4 = boost::make_unique<BDABuffer>(kDataSize);
  auto buffer5 = boost::make_unique<BDABuffer>(kDataSize);
  AddShortRow(*buffer1, kTime + 0.0, 41, false);
  AddShortRow(*buffer2, kTime + 10.0, 42, false);
  AddShortRow(*buffer3, kTime + 20.0, 43, false);
  AddShortRow(*buffer4, kTime + 30.0, 44, false);
  AddShortRow(*buffer5, kTime + 40.0, 45, false);

  // An interval is complete if the interval buffer has a row that starts
  // at or after <end_of_interval> + <max_row_interval>.

  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer1));
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer2));
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer3));
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer4));
  BOOST_CHECK(interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 4u);
  CheckBuffer(interval, *buffer1);

  interval.Advance(kInterval);
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 3u);

  BOOST_CHECK_NO_THROW(interval.AddBuffer(*buffer5));
  BOOST_CHECK(interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 4u);
  CheckBuffer(interval, *buffer2);

  interval.Advance(kInterval);
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 3u);
  CheckBuffer(interval, *buffer3);

  interval.Advance(kInterval);
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 2u);
  CheckBuffer(interval, *buffer4);

  interval.Advance(kInterval);
  BOOST_CHECK(!interval.IsComplete());
  BOOST_CHECK_EQUAL(interval.GetPendingBufferCount(), 1u);
  CheckBuffer(interval, *buffer5);

  CheckAdvanceBeyondEnd(interval);
}

BOOST_AUTO_TEST_SUITE_END()
