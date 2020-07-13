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
/// @brief Unit tests for the EnumBitset template class.
/// @author Maik Nijhuis

#include "../../EnumBitset.h"

#include <boost/test/unit_test.hpp>

using DP3::EnumBitset;

namespace
{
  enum class TestEnum { kBit0, kBit1, kBit2, kBitCount };
}

BOOST_AUTO_TEST_SUITE( enumbitset )

BOOST_AUTO_TEST_CASE( default_constructor )
{
  EnumBitset<TestEnum> bitset;
  BOOST_CHECK(!bitset[TestEnum::kBit0]);
  BOOST_CHECK(!bitset[TestEnum::kBit1]);
  BOOST_CHECK(!bitset[TestEnum::kBit2]);
}

BOOST_AUTO_TEST_CASE( list_constructor )
{
  EnumBitset<TestEnum> bitset{true, false, true};
  BOOST_CHECK(bitset[TestEnum::kBit0]);
  BOOST_CHECK(!bitset[TestEnum::kBit1]);
  BOOST_CHECK(bitset[TestEnum::kBit2]);
}

BOOST_AUTO_TEST_CASE( set_single_bit )
{
  EnumBitset<TestEnum> bitset;
  BOOST_CHECK(!bitset[TestEnum::kBit1]);
  bitset[TestEnum::kBit1] = true;
  BOOST_CHECK(bitset[TestEnum::kBit1]);
}

BOOST_AUTO_TEST_CASE( set_all_bits )
{
  EnumBitset<TestEnum> bitset;
  BOOST_CHECK_EQUAL(&bitset.Set(), &bitset);
  BOOST_CHECK(bitset[TestEnum::kBit0]);
  BOOST_CHECK(bitset[TestEnum::kBit1]);
  BOOST_CHECK(bitset[TestEnum::kBit2]);
}

BOOST_AUTO_TEST_SUITE_END()
