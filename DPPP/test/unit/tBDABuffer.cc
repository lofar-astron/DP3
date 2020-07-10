// tScaleData.cc: Test program for class BDABuffer
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
//
// @author Lars Krombeen

#include "../Common/Types.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../../BDABuffer.h"

using DP3::DPPP::BDABuffer;

BOOST_AUTO_TEST_SUITE(bdabuffer)

BOOST_AUTO_TEST_CASE( initialization )
{
    BDABuffer buffer {2};
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {0});
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
}

BOOST_AUTO_TEST_CASE( copy )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};
    BDABuffer buffer {4 * n_channels * n_correlations};
    const bool flags[n_channels * n_correlations] {true};
    const std::complex<float> data {1};
    const float weights[n_channels * n_correlations] {1};

    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data);
    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags);
    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights);

    BDABuffer buffer_copy {buffer};

    BOOST_CHECK(buffer.GetData() != buffer_copy.GetData());
    BOOST_CHECK(buffer.GetData() != nullptr);
    BOOST_CHECK(buffer.GetFlags() != nullptr);
    BOOST_CHECK(buffer.GetWeights() != nullptr);
    BOOST_CHECK(buffer.GetFullResFlags() != nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetData(), *buffer_copy.GetData());
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), *buffer_copy.GetFlags());
    BOOST_CHECK(std::isnan(*buffer.GetWeights()) && std::isnan(*buffer_copy.GetWeights()));
    BOOST_CHECK_EQUAL(*buffer.GetFullResFlags(), *buffer_copy.GetFullResFlags());
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), buffer_copy.GetNumberOfElements());
    BOOST_CHECK(buffer_copy.AddRow(1., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(!buffer_copy.AddRow(1., 1., 0, 0, n_channels, n_correlations, &data));
}

BOOST_AUTO_TEST_CASE( wrong_add_row_order )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};
    BDABuffer buffer {3 * n_channels * n_correlations};
    const std::complex<float> data {1};

    buffer.AddRow(1., 1., 0, 0, n_channels, n_correlations, &data);
    BOOST_CHECK_THROW(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( add_null_data )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, true, false, false, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK(buffer.GetData() != nullptr);
    BOOST_CHECK(std::isnan((*buffer.GetData()).real()));
    BOOST_CHECK(std::isnan((*buffer.GetData()).imag()));
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
}

BOOST_AUTO_TEST_CASE( add_null_flag )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, false, true, false, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), false);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
}

BOOST_AUTO_TEST_CASE( add_null_weight )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, false, false, true, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK(buffer.GetWeights() != nullptr);
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
    BOOST_CHECK(std::isnan(*buffer.GetWeights()));
}

BOOST_AUTO_TEST_CASE( data )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations};
    const std::complex<float> data {1};
    const std::complex<float> data2 {2};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data2));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {n_channels * n_correlations});
    BOOST_CHECK(buffer.GetData() != nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetData(), data);
    BOOST_CHECK_EQUAL(*(buffer.GetData() + 1), data2);
}

BOOST_AUTO_TEST_CASE( flags )
{
    const size_t n_channels {1};
    const size_t n_correlations {5};
    BDABuffer buffer {n_channels * n_correlations};
    const bool flags[n_channels * n_correlations] {true};
    const bool flags2[n_channels * n_correlations] {false};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags2));
    BOOST_CHECK(buffer.GetFlags() != nullptr);
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {n_channels * n_correlations});
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), flags[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetFlags() + 1), flags2[0]);
}

BOOST_AUTO_TEST_CASE( weights )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    // Add 1 to check if a 3rd does not get added if there is 1 pool left
    BDABuffer buffer {2 * n_channels * n_correlations + 1};
    const float weights[] {1, 2, 3, 4, 5, 6};
    const float weights2[] {3, 4 ,5 ,6 ,7 ,8};
    const float weights3[] {5, 6, 7, 8, 9, 10};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights2));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights3));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {2 * n_channels * n_correlations});
    BOOST_CHECK(buffer.GetWeights() != nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetWeights(), weights[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 1), weights[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 2), weights2[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 3), weights2[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 4), weights3[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 5), weights3[1]);
}

BOOST_AUTO_TEST_CASE( full_res_flags )
{
    const size_t n_channels {1};
    const size_t n_correlations {5};
    BDABuffer buffer {n_channels * n_correlations};
    const bool flags[n_channels * n_correlations] {true};
    const bool flags2[n_channels * n_correlations] {false};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, nullptr, flags));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, nullptr, flags2));
    BOOST_CHECK(buffer.GetFullResFlags() != nullptr);
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {n_channels * n_correlations});
    BOOST_CHECK_EQUAL(*buffer.GetFullResFlags(), flags[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetFullResFlags() + 1), flags2[0]);
}

BOOST_AUTO_TEST_CASE( uvw )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    BDABuffer buffer {2 * n_channels * n_correlations};
    const double uvw[n_channels * n_correlations] {1, 2, 3};
    const double uvw2[n_channels * n_correlations] {3, 4, 5};
    const double uvw3[n_channels * n_correlations] {5, 6, 7};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, nullptr, nullptr, uvw));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, nullptr, nullptr, uvw2));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, nullptr, nullptr, uvw3));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {2 * n_channels * n_correlations});
    const std::vector<BDABuffer::Row>& rows = buffer.GetRows();
    BOOST_CHECK_EQUAL(rows[0].uvw_[0], uvw[0]);
    BOOST_CHECK_EQUAL(rows[0].uvw_[1], uvw[1]);
    BOOST_CHECK_EQUAL(rows[1].uvw_[0], uvw2[0]);
    BOOST_CHECK_EQUAL(rows[1].uvw_[1], uvw2[1]);
}

BOOST_AUTO_TEST_CASE( check_row_data_size )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    BDABuffer buffer {n_channels * n_correlations};
    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations);

    const std::vector<BDABuffer::Row>& rows = buffer.GetRows();
    BOOST_CHECK_EQUAL(rows[0].GetDataSize(), n_channels * n_correlations);
}

BOOST_AUTO_TEST_CASE( clear )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    // Add 1 to check if a 4th does not get added if there is 1 pool left
    BDABuffer buffer {3 * n_channels * n_correlations + 1};
    const bool flags[n_channels * n_correlations] {true};
    const bool fr_flags[n_channels * n_correlations] {true, true};
    const std::complex<float> data {1};
    const float weights[n_channels * n_correlations] {1};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights, fr_flags));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {3 * n_channels * n_correlations});
    buffer.Clear();
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {0});

    // Check that 3 elements can be added again
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
}

BOOST_AUTO_TEST_CASE( get_rows )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    // Add 1 to check if a 4th does not get added if there is 1 pool left
    BDABuffer buffer {3 * n_channels * n_correlations + 1};
    const bool flags[n_channels * n_correlations] {true};
    const bool fr_flags[n_channels * n_correlations] {false};
    const std::complex<float> data {1};
    const float weights[n_channels * n_correlations] {1};

    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data, flags, weights, fr_flags);

    BOOST_CHECK(buffer.GetData() != nullptr);
    BOOST_CHECK(buffer.GetFlags() != nullptr);
    BOOST_CHECK(buffer.GetWeights() != nullptr);
    BOOST_CHECK(buffer.GetFullResFlags() != nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetData(0), data);
    BOOST_CHECK_EQUAL(*buffer.GetFlags(0), flags[0]);
    BOOST_CHECK_EQUAL(*buffer.GetFullResFlags(0), fr_flags[0]);
    BOOST_CHECK_EQUAL(*buffer.GetWeights(0), weights[0]);
}

BOOST_AUTO_TEST_CASE( time_is_less ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(5., 10.), true);
}

BOOST_AUTO_TEST_CASE( time_is_less_false ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(15., 10.), false);
}

BOOST_AUTO_TEST_CASE( time_is_less_near ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(9.99999998, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( time_is_less_near_false ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsLess(9.99999999, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( time_is_greater_false ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(5., 10.), false);
}

BOOST_AUTO_TEST_CASE( time_is_greater ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(15., 10.), true);
}

BOOST_AUTO_TEST_CASE( time_is_greater_near_false ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(9.99999998, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( time_is_greater_near ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(9.999999999, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( time_is_greater_equal ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsGreaterEqual(10.0000000, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( time_is_equal ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(10.0000000, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( time_is_not_equal ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(20.0000000, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( time_is_near_equal ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(9.999999991, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( time_is_not_near_equal ) {
    BOOST_CHECK_EQUAL(BDABuffer::TimeIsEqual(9.99999999, 10.0000000), false);
}

BOOST_AUTO_TEST_SUITE_END()
