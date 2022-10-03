// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Fields.h"

#include <sstream>

#include <boost/test/unit_test.hpp>

using dp3::common::Fields;

BOOST_AUTO_TEST_SUITE(fields)

namespace {
void CheckFields(const Fields& fields, bool data, bool flags, bool weights,
                 bool full_res_flags, bool uvw) {
  BOOST_CHECK(fields.Data() == data);
  BOOST_CHECK(fields.Flags() == flags);
  BOOST_CHECK(fields.Weights() == weights);
  BOOST_CHECK(fields.FullResFlags() == full_res_flags);
  BOOST_CHECK(fields.Uvw() == uvw);
}
}  // namespace

BOOST_AUTO_TEST_CASE(constructor_default) {
  const Fields fields;
  CheckFields(fields, false, false, false, false, false);
}

BOOST_AUTO_TEST_CASE(constructor_single_field) {
  const Fields data(Fields::Single::kData);
  const Fields flags(Fields::Single::kFlags);
  const Fields weights(Fields::Single::kWeights);
  const Fields full_res_flags(Fields::Single::kFullResFlags);
  const Fields uvw(Fields::Single::kUvw);
  CheckFields(data, true, false, false, false, false);
  CheckFields(flags, false, true, false, false, false);
  CheckFields(weights, false, false, true, false, false);
  CheckFields(full_res_flags, false, false, false, true, false);
  CheckFields(uvw, false, false, false, false, true);
}

BOOST_AUTO_TEST_CASE(combine) {
  Fields fields;
  fields |= Fields(Fields::Single::kFlags);
  CheckFields(fields, false, true, false, false, false);
  fields |= Fields(Fields::Single::kFullResFlags);
  CheckFields(fields, false, true, false, true, false);

  const Fields data(Fields::Single::kData);
  const Fields uvw(Fields::Single::kUvw);
  CheckFields(data | uvw, true, false, false, false, true);
}

BOOST_AUTO_TEST_CASE(equality_operators) {
  const Fields left_empty;
  const Fields right_empty;
  const Fields left_three = Fields(Fields::Single::kData) |
                            Fields(Fields::Single::kFlags) |
                            Fields(Fields::Single::kUvw);
  const Fields right_three = Fields(Fields::Single::kData) |
                             Fields(Fields::Single::kFlags) |
                             Fields(Fields::Single::kUvw);
  // Check operator==
  BOOST_CHECK(left_empty == left_empty);
  BOOST_CHECK(left_empty == right_empty);
  BOOST_CHECK(right_three == right_three);
  BOOST_CHECK(left_three == right_three);
  BOOST_CHECK(!(left_empty == right_three));
  BOOST_CHECK(!(right_three == left_empty));

  // Check operator!=
  BOOST_CHECK(!(left_empty != left_empty));
  BOOST_CHECK(!(left_empty != right_empty));
  BOOST_CHECK(!(right_three != right_three));
  BOOST_CHECK(!(left_three != right_three));
  BOOST_CHECK(left_empty != right_three);
  BOOST_CHECK(right_three != left_empty);
}

BOOST_AUTO_TEST_CASE(outputstream) {
  std::stringstream stream_empty;
  stream_empty << Fields();
  BOOST_CHECK_EQUAL(stream_empty.str(), "[]");

  std::stringstream stream_single;
  stream_single << Fields(Fields::Single::kData);
  BOOST_CHECK_EQUAL(stream_single.str(), "[data]");

  std::stringstream stream_all;
  stream_all << (Fields(Fields::Single::kData) |
                 Fields(Fields::Single::kWeights) |
                 Fields(Fields::Single::kFlags) |
                 Fields(Fields::Single::kFullResFlags) |
                 Fields(Fields::Single::kUvw));
  BOOST_CHECK_EQUAL(stream_all.str(), "[data, flags, weights, frf, uvw]");
}

BOOST_AUTO_TEST_SUITE_END()
