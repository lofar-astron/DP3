// tScaleData.cc: Test program for class ScaleData
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

BOOST_AUTO_TEST_SUITE(dbabuffer)

BOOST_AUTO_TEST_CASE( test_bda_initialization )
{
    BDABuffer buffer {2};
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {0});
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFullResFlags(), nullptr);
}

BOOST_AUTO_TEST_CASE( test_bda_add_null_data )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, true, false, false, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK(std::isnan((*buffer.GetData()).real()));
    BOOST_CHECK(std::isnan((*buffer.GetData()).imag()));
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
}

BOOST_AUTO_TEST_CASE( test_bda_add_null_flag )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, false, true, false, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), false);
    BOOST_CHECK_EQUAL(buffer.GetWeights(), nullptr);
}

BOOST_AUTO_TEST_CASE( test_bda_add_null_weight )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations, false, false, true, false};
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations));
    BOOST_CHECK_EQUAL(buffer.GetData(), nullptr);
    BOOST_CHECK_EQUAL(buffer.GetFlags(), nullptr);
    BOOST_CHECK(std::isnan(*buffer.GetWeights()));
}

BOOST_AUTO_TEST_CASE( test_dba_data )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations};
    const std::complex<float> data {1};
    const std::complex<float> data2 {2};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data2));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {n_channels * n_correlations});
    BOOST_CHECK_EQUAL(*buffer.GetData(), data);
    BOOST_CHECK_EQUAL(*(buffer.GetData() + 1), data2);
}

BOOST_AUTO_TEST_CASE( test_dba_flags )
{
    const size_t n_channels {1};
    const size_t n_correlations {5};
    BDABuffer buffer {n_channels * n_correlations};
    bool flags[] {true};
    bool flags2[] {false};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags2));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {n_channels * n_correlations});
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), flags[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetFlags() + 1), flags2[0]);
}

// TODO test that memory allocation is correct in constructor

BOOST_AUTO_TEST_CASE( test_dba_weights )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};

    BDABuffer buffer {2 * n_channels * n_correlations + 1};
    float weights[] {1, 2};
    float weights2[] {3, 4};
    float weights3[] {5, 6};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights2));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights3));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {2 * n_channels * n_correlations});
    BOOST_CHECK_EQUAL(*buffer.GetWeights(), weights[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 1), weights[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 2), weights2[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 3), weights2[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 4), weights3[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 5), weights3[1]);
}

BOOST_AUTO_TEST_CASE( test_dba_clear )
{
    const size_t n_channels {2};
    const size_t n_correlations {3};
    BDABuffer buffer {3 * n_channels * n_correlations + 1};
    bool flags[] {true};
    const std::complex<float> data {1};
    float weights[] {1};

    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags));
    BOOST_CHECK(buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights));
    BOOST_CHECK(!buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data));
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {3 * n_channels * n_correlations});
    buffer.Clear();
    BOOST_CHECK_EQUAL(buffer.GetNumberOfElements(), size_t {0});
}

BOOST_AUTO_TEST_CASE( test_modify_data ) 
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations};
    const std::complex<float> data {1};

    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, &data);

    const std::complex<float> new_data {2};
    *buffer.GetData() = new_data;
    BOOST_CHECK_EQUAL(*buffer.GetData(), new_data);
}

BOOST_AUTO_TEST_CASE( test_modify_flags )
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations};
    bool flags[] {false, true, false};

    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, flags);
    
    bool new_flags[] {true, false, true};
    *buffer.GetFlags() = new_flags[0];
    *(buffer.GetFlags() + 1) = new_flags[1];
    *(buffer.GetFlags() + 2) = new_flags[2];
    BOOST_CHECK_EQUAL(*buffer.GetFlags(), new_flags[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetFlags() + 1), new_flags[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetFlags() + 2), new_flags[2]);
}

BOOST_AUTO_TEST_CASE( test_modify_weights ) 
{
    const size_t n_channels {2};
    const size_t n_correlations {4};
    BDABuffer buffer {n_channels * n_correlations};
    float weights[] {1, 2, 3};

    buffer.AddRow(0., 1., 0, 0, n_channels, n_correlations, nullptr, nullptr, weights);

    float new_weights[] {4, 5, 6};
    *buffer.GetWeights() = new_weights[0];
    *(buffer.GetWeights() + 1) = new_weights[1];
    *(buffer.GetWeights() + 2) = new_weights[2];
    BOOST_CHECK_EQUAL(*buffer.GetWeights(), new_weights[0]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 1), new_weights[1]);
    BOOST_CHECK_EQUAL(*(buffer.GetWeights() + 2), new_weights[2]);
}

BOOST_AUTO_TEST_CASE( test_time_is_less ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsLess(5., 10.), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_less_false ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsLess(15., 10.), false);
}

BOOST_AUTO_TEST_CASE( test_time_is_less_near ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsLess(9.99999998, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_less_near_false ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsLess(9.99999999, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( test_time_is_greater_false ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsGreaterEqual(5., 10.), false);
}

BOOST_AUTO_TEST_CASE( test_time_is_greater ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsGreaterEqual(15., 10.), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_greater_near_false ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsGreaterEqual(9.99999998, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( test_time_is_greater_near ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsGreaterEqual(9.999999999, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_greater_equal ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsGreaterEqual(10.0000000, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_equal ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsEqual(10.0000000, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_not_equal ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsEqual(20.0000000, 10.0000000), false);
}

BOOST_AUTO_TEST_CASE( test_time_is_near_equal ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsEqual(9.999999991, 10.0000000), true);
}

BOOST_AUTO_TEST_CASE( test_time_is_not_near_equal ) {
    BDABuffer buffer {0};
    BOOST_CHECK_EQUAL(buffer.TimeIsEqual(9.99999999, 10.0000000), false);
}

// TODO add tests for the different flags  if data is reserved or not

BOOST_AUTO_TEST_SUITE_END()
