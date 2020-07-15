// Copyright (C) 2020
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
/// @brief Unit tests for the BDABuffer class.
/// @author Lars Krombeen & Maik Nijhuis

#include "../Common/Types.h"

#include <boost/test/unit_test.hpp>

#include "../../BDABuffer.h"

using DP3::DPPP::BDABuffer;

namespace {
    const double kTime = {0.0};
    const double kInterval {1.0};
    const DP3::rownr_t kRowNr {1};
    const std::size_t kBaselineNr {42};
    const std::size_t kNChannels {2};
    const std::size_t kNCorrelations {3};
    const std::size_t kDataSize {kNChannels * kNCorrelations};
    const std::size_t k1Channel {1};
    const std::size_t k1Correlation {1};
    const std::size_t k1DataSize {k1Channel * k1Correlation};

    const double kTimeEpsilon = 1.0e-8;
    const double kHalfEpsilon = kTimeEpsilon / 2;
    const double kTwoEpsilon = 2.0 * kTimeEpsilon;
}

BOOST_AUTO_TEST_SUITE(bdabuffer)

BOOST_AUTO_TEST_CASE( initialization )
{
    BDABuffer buffer {2};
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 0u);
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetRows().size(), 0u);
}

BOOST_AUTO_TEST_CASE( copy )
{
    const std::complex<float> kData1[kDataSize] {
        {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}
    };
    const std::complex<float> kData2[kDataSize] {
        {-1, -1}, {-2, -2}, {-3 ,-3}, {-4, -4}, {-5, -5}, {-6, -6}
    };
    const bool kFlags[kDataSize] {true, false, true, true, false ,true};
    const float kWeights1[kDataSize] {21, 22, 23, 24, 25, 26};
    const float kWeights2[kDataSize] {31, 32, 33, 34, 35, 36};
    const double kUvw[3] {41, 42, 43};

    BDABuffer buffer {3 * kDataSize};
    buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr, kNChannels,
                  kNCorrelations, kData1, nullptr, kWeights1, kFlags, kUvw);
    buffer.AddRow(kTime + 1., kInterval + 1., kRowNr + 1, kBaselineNr + 1, kNChannels,
                  kNCorrelations, kData2, kFlags, kWeights2, nullptr, kUvw);

    BDABuffer buffer_copy {buffer};

    // Verify the memory pool data in the copy.
    BOOST_CHECK(buffer.GetData() != buffer_copy.GetData());
    BOOST_CHECK(buffer.GetFlags() != buffer_copy.GetFlags());
    BOOST_CHECK(buffer.GetWeights() != buffer_copy.GetWeights());
    BOOST_CHECK(buffer.GetFullResFlags() != buffer_copy.GetFullResFlags());
    BOOST_CHECK(buffer.GetData() != nullptr);
    BOOST_CHECK(buffer.GetFlags() != nullptr);
    BOOST_CHECK(buffer.GetWeights() != nullptr);
    BOOST_CHECK(buffer.GetFullResFlags() != nullptr);

    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), buffer_copy.GetNumberOfElements());
    for (std::size_t i = 0; i < buffer.GetNumberOfElements(); ++i) {
        BOOST_CHECK_EQUAL(buffer.GetData()[i], buffer_copy.GetData()[i]);
        BOOST_CHECK_EQUAL(buffer.GetFlags()[i], buffer_copy.GetFlags()[i]);
        BOOST_CHECK_EQUAL(buffer.GetWeights()[i], buffer_copy.GetWeights()[i]);
        BOOST_CHECK_EQUAL(buffer.GetFullResFlags()[i], buffer_copy.GetFullResFlags()[i]);
    }

    // Verify the copied rows.
    const auto& rows_copy = buffer_copy.GetRows();
    BOOST_CHECK_EQUAL(rows_copy.size(), 2u);
    BOOST_CHECK_EQUAL(rows_copy[0].GetDataSize(), kDataSize);
    BOOST_CHECK_EQUAL(rows_copy[0].time_, kTime);
    BOOST_CHECK_EQUAL(rows_copy[0].interval_, kInterval);
    BOOST_CHECK_EQUAL(rows_copy[0].row_nr_, kRowNr);
    BOOST_CHECK_EQUAL(rows_copy[0].baseline_nr_, kBaselineNr);
    BOOST_CHECK_EQUAL(rows_copy[0].n_channels_, kNChannels);
    BOOST_CHECK_EQUAL(rows_copy[0].n_correlations_, kNCorrelations);

    BOOST_CHECK_EQUAL(rows_copy[1].GetDataSize(), kDataSize);
    BOOST_CHECK_EQUAL(rows_copy[1].time_, kTime + 1.);
    BOOST_CHECK_EQUAL(rows_copy[1].interval_, kInterval + 1.);
    BOOST_CHECK_EQUAL(rows_copy[1].row_nr_, kRowNr + 1);
    BOOST_CHECK_EQUAL(rows_copy[1].baseline_nr_, kBaselineNr + 1);
    BOOST_CHECK_EQUAL(rows_copy[1].n_channels_, kNChannels);
    BOOST_CHECK_EQUAL(rows_copy[1].n_correlations_, kNCorrelations);

    for (std::size_t i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(rows_copy[0].uvw_[i], kUvw[i]);
        BOOST_CHECK_EQUAL(rows_copy[1].uvw_[i], kUvw[i]);
    }

    // Verify that the original has remaining capacity, but the copy doesn't.
    BOOST_CHECK(buffer.AddRow(kTime + 2., kInterval, kRowNr, kBaselineNr + 2,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(!buffer_copy.AddRow(kTime + 2., kInterval, kRowNr, kBaselineNr + 2,
                                    kNChannels, kNCorrelations));
}

BOOST_AUTO_TEST_CASE( add_all_fields )
{
    const std::complex<float> kData1[kDataSize] {
        {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}
    };
    const std::complex<float> kData2[1] { {7, 7} };
    const bool kFlags1[kDataSize] {true, false, true, true, false, true};
    const bool kFlags2[k1DataSize] {false };
    const float kWeights1[kDataSize] {21, 22, 23, 24, 25, 26};
    const float kWeights2[k1DataSize] {27};
    const double kUvw1[3] {31, 32, 33};
    const double kUvw2[3] {41, 42, 43};

    BDABuffer buffer {kDataSize + k1DataSize};

    // Add rows and verify that the capacity is reached at the third row.
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations,
                              kData1, kFlags1, kWeights1, kFlags1, kUvw1));
    BOOST_CHECK(buffer.AddRow(kTime + 1., kInterval + 1., kRowNr + 1, kBaselineNr + 1,
                              k1Channel, k1Correlation,
                              kData2, kFlags2, kWeights2, kFlags2, kUvw2));
    BOOST_CHECK(!buffer.AddRow(kTime + 2., kInterval, kRowNr, kBaselineNr + 2,
                               k1Channel, k1Correlation));

    // Verify the data in the buffer.
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), kDataSize + k1DataSize);
    BOOST_CHECK(buffer.GetData());
    BOOST_CHECK(buffer.GetFlags());
    BOOST_CHECK(buffer.GetWeights());
    BOOST_CHECK(buffer.GetFullResFlags());

    BOOST_CHECK_EQUAL(buffer.GetData(0), buffer.GetData());
    BOOST_CHECK_EQUAL(buffer.GetData(1), buffer.GetData() + kDataSize);
    BOOST_CHECK_EQUAL(buffer.GetFlags(0), buffer.GetFlags());
    BOOST_CHECK_EQUAL(buffer.GetFlags(1), buffer.GetFlags() + kDataSize);
    BOOST_CHECK_EQUAL(buffer.GetWeights(0), buffer.GetWeights());
    BOOST_CHECK_EQUAL(buffer.GetWeights(1), buffer.GetWeights() + kDataSize);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(0), buffer.GetFullResFlags());
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(1), buffer.GetFullResFlags() + kDataSize);

    for (std::size_t i = 0; i < kDataSize; ++i) {
        BOOST_CHECK_EQUAL(buffer.GetData()[i], kData1[i]);
        BOOST_CHECK_EQUAL(buffer.GetFlags()[i], kFlags1[i]);
        BOOST_CHECK_EQUAL(buffer.GetWeights()[i], kWeights1[i]);
        BOOST_CHECK_EQUAL(buffer.GetFullResFlags()[i], kFlags1[i]);
    }
    for (std::size_t i = 0; i < k1DataSize; ++i) {
        BOOST_CHECK_EQUAL(buffer.GetData(1)[i], kData2[i]);
        BOOST_CHECK_EQUAL(buffer.GetFlags(1)[i], kFlags2[i]);
        BOOST_CHECK_EQUAL(buffer.GetWeights(1)[i], kWeights2[i]);
        BOOST_CHECK_EQUAL(buffer.GetFullResFlags(1)[i], kFlags2[i]);
    }

    // Verify the rows.
    const auto& rows = buffer.GetRows();
    BOOST_CHECK_EQUAL(rows.size(), 2u);
    BOOST_CHECK_EQUAL(rows[0].GetDataSize(), kDataSize);
    BOOST_CHECK_EQUAL(rows[0].time_, kTime);
    BOOST_CHECK_EQUAL(rows[0].interval_, kInterval);
    BOOST_CHECK_EQUAL(rows[0].row_nr_, kRowNr);
    BOOST_CHECK_EQUAL(rows[0].baseline_nr_, kBaselineNr);
    BOOST_CHECK_EQUAL(rows[0].n_channels_, kNChannels);
    BOOST_CHECK_EQUAL(rows[0].n_correlations_, kNCorrelations);

    BOOST_CHECK_EQUAL(rows[1].GetDataSize(), k1DataSize);
    BOOST_CHECK_EQUAL(rows[1].time_, kTime + 1.);
    BOOST_CHECK_EQUAL(rows[1].interval_, kInterval + 1.);
    BOOST_CHECK_EQUAL(rows[1].row_nr_, kRowNr + 1);
    BOOST_CHECK_EQUAL(rows[1].baseline_nr_, kBaselineNr + 1);
    BOOST_CHECK_EQUAL(rows[1].n_channels_, k1Channel);
    BOOST_CHECK_EQUAL(rows[1].n_correlations_, k1Correlation);

    for (std::size_t i = 0; i < 2; ++i) {
        BOOST_CHECK_EQUAL(rows[i].data_, buffer.GetData(i));
        BOOST_CHECK_EQUAL(rows[i].flags_, buffer.GetFlags(i));
        BOOST_CHECK_EQUAL(rows[i].weights_, buffer.GetWeights(i));
        BOOST_CHECK_EQUAL(rows[i].full_res_flags_, buffer.GetFullResFlags(i));
    }

    for (std::size_t i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(rows[0].uvw_[i], kUvw1[i]);
        BOOST_CHECK_EQUAL(rows[1].uvw_[i], kUvw2[i]);
    }
}

BOOST_AUTO_TEST_CASE( add_no_fields )
{
    BDABuffer buffer {kDataSize};
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(buffer.GetData());
    BOOST_CHECK(buffer.GetFlags());
    BOOST_CHECK(buffer.GetWeights());
    BOOST_CHECK(buffer.GetFullResFlags());

    // Check that the memory pools hold default values.
    for (std::size_t i = 0; i < kDataSize; ++i) {
        BOOST_CHECK(std::isnan(buffer.GetData()[i].real()));
        BOOST_CHECK(std::isnan(buffer.GetData()[i].imag()));
        BOOST_CHECK_EQUAL(buffer.GetFlags()[i], false);
        BOOST_CHECK(std::isnan(buffer.GetWeights()[i]));
        BOOST_CHECK_EQUAL(buffer.GetFullResFlags()[i], false);
    }

    // Verify the row.
    const auto& rows = buffer.GetRows();
    BOOST_CHECK_EQUAL(rows.size(), 1u);
    const auto& row = rows.front();
    BOOST_CHECK_EQUAL(row.GetDataSize(), kDataSize);
    BOOST_CHECK_EQUAL(row.time_, kTime);
    BOOST_CHECK_EQUAL(row.interval_, kInterval);
    BOOST_CHECK_EQUAL(row.row_nr_, kRowNr);
    BOOST_CHECK_EQUAL(row.baseline_nr_, kBaselineNr);
    BOOST_CHECK_EQUAL(row.n_channels_, kNChannels);
    BOOST_CHECK_EQUAL(row.n_correlations_, kNCorrelations);
    BOOST_CHECK_EQUAL(row.data_, buffer.GetData(0));
    BOOST_CHECK_EQUAL(row.flags_, buffer.GetFlags(0));
    BOOST_CHECK_EQUAL(row.weights_, buffer.GetWeights(0));
    BOOST_CHECK_EQUAL(row.full_res_flags_, buffer.GetFullResFlags(0));
    for (std::size_t i = 0; i < 3; ++i) {
        BOOST_CHECK(std::isnan(row.uvw_[i]));
    }
}

BOOST_AUTO_TEST_CASE( disabled_fields )
{
    const std::complex<float> kData[kDataSize] { {1, 2} };
    const bool kFlags[kDataSize] {true};
    const float kWeights[kDataSize] {42};
    const double kUvw[3] {42};

    BDABuffer::Fields fields;
    fields.data_ = false;
    fields.flags_ = false;
    fields.weights_ = false;
    fields.full_res_flags_ = false;

    BDABuffer buffer {kDataSize, fields};
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations,
                              kData, kFlags, kWeights, kFlags, kUvw));

    // Verify the data pointers from buffer.Get*().
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);

    BOOST_CHECK_EQUAL(buffer.GetData(0), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(0), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(0), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(0), nullptr);

    // Verify the data pointers in the row.
    const auto& rows = buffer.GetRows();
    BOOST_CHECK_EQUAL(rows.size(), 1u);
    const auto& row = rows.front();
    BOOST_CHECK_EQUAL(row.data_, nullptr);
    BOOST_CHECK_EQUAL(row.flags_, nullptr);
    BOOST_CHECK_EQUAL(row.weights_, nullptr);
    BOOST_CHECK_EQUAL(row.full_res_flags_, nullptr);
}

BOOST_AUTO_TEST_CASE( add_wrong_ordering )
{
    const double kTime1 = 3.0;
    const double kInterval1 = 3.0;
    const double kTime2 = 0.0;
    const double kInterval2 = 3.0;

    BDABuffer buffer {2 * kDataSize};

    BOOST_CHECK(buffer.AddRow(kTime1, kInterval1, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations));
    // end time equals start time of last row -> invalid ordering.
    BOOST_CHECK_THROW(buffer.AddRow(kTime2, kInterval2, kRowNr, kBaselineNr,
                                    kNChannels, kNCorrelations),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( add_overlap )
{
    const double kTime1 = 3.0;
    const double kInterval1 = 3.0;
    const double kTimePartialOverlap = 2.0;
    const double kIntervalPartialOverlap = 2.0;
    const double kTimeFullOverlap = 1.0;
    const double kIntervalFullOverlap = 10.0;

    BDABuffer buffer {3 * kDataSize};

    BOOST_CHECK(buffer.AddRow(kTime1, kInterval1, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations));
    // Add partially overlapping row, starting before the last row.
    BOOST_CHECK(buffer.AddRow(kTimePartialOverlap, kIntervalPartialOverlap,
                              kRowNr, kBaselineNr, kNChannels, kNCorrelations));
    // Add fully overlapping row, starting before the last row.
    BOOST_CHECK(buffer.AddRow(kTimeFullOverlap, kIntervalFullOverlap,
                              kRowNr, kBaselineNr, kNChannels, kNCorrelations));
}

BOOST_AUTO_TEST_CASE( clear )
{
    BDABuffer buffer {3 * kDataSize};

    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 1,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 2,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(!buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 3,
                               kNChannels, kNCorrelations));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 3 * kDataSize);
    BOOST_CHECK_EQUAL(buffer.GetRows().size(), 3u);

    buffer.Clear();
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), 0u);
    BOOST_CHECK_EQUAL(buffer.GetRows().size(), 0u);

    // Check that 3 rows can be added again.
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 1,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 2,
                              kNChannels, kNCorrelations));
    BOOST_CHECK(!buffer.AddRow(kTime, kInterval, kRowNr, kBaselineNr + 3,
                               kNChannels, kNCorrelations));
}

BOOST_AUTO_TEST_CASE( time_is_less )
{
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(5., 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(15., 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(10. - kTwoEpsilon, 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(10. - kHalfEpsilon, 10.), false);
}

BOOST_AUTO_TEST_CASE( time_is_less_equal )
{
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLessEqual(5., 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLessEqual(15., 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLessEqual(10. + kTwoEpsilon, 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLessEqual(10. + kHalfEpsilon, 10.), true);
}

BOOST_AUTO_TEST_CASE( time_is_greater_equal )
{
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(5., 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(15., 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(10., 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(10. - kTwoEpsilon, 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(10. - kHalfEpsilon, 10.), true);
}

BOOST_AUTO_TEST_CASE( time_is_equal )
{
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(10., 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(20., 10.), false);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(10. + kHalfEpsilon, 10.), true);
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(10. - kTwoEpsilon, 10.), false);
}

BOOST_AUTO_TEST_SUITE_END()
