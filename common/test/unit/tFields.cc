// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <dp3/common/Fields.h>

#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>

using dp3::common::Fields;

BOOST_AUTO_TEST_SUITE(fields)

namespace {

void CheckFields(const Fields& fields, std::string set_fields) {
  BOOST_CHECK(fields.Data() == (set_fields.find("data") != std::string::npos));
  BOOST_CHECK(fields.Flags() ==
              (set_fields.find("flags") != std::string::npos));
  BOOST_CHECK(fields.Weights() ==
              (set_fields.find("weights") != std::string::npos));
  BOOST_CHECK(fields.FullResFlags() ==
              (set_fields.find("fullResFlags") != std::string::npos));
  BOOST_CHECK(fields.Uvw() == (set_fields.find("uvw") != std::string::npos));
}

}  // namespace

BOOST_AUTO_TEST_CASE(constructor_default) {
  const Fields fields;
  CheckFields(fields, "");
}

BOOST_AUTO_TEST_CASE(constructor_single_field) {
  const Fields data(Fields::Single::kData);
  const Fields flags(Fields::Single::kFlags);
  const Fields weights(Fields::Single::kWeights);
  const Fields full_res_flags(Fields::Single::kFullResFlags);
  const Fields uvw(Fields::Single::kUvw);
  CheckFields(data, "data");
  CheckFields(flags, "flags");
  CheckFields(weights, "weights");
  CheckFields(full_res_flags, "fullResFlags");
  CheckFields(uvw, "uvw");
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
  BOOST_CHECK_EQUAL(stream_all.str(),
                    "[data, flags, weights, fullresflags, uvw]");
}

BOOST_AUTO_TEST_CASE(update_requirements) {
  const Fields kDataField(Fields::Single::kData);

  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(), Fields(Fields::Single::kData)),
              "");
  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(Fields::Single::kData),
                                      Fields(Fields::Single::kData)),
              "data");
  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(Fields::Single::kData), Fields()),
              "data");
  CheckFields(Fields(kDataField).UpdateRequirements(Fields(), Fields()),
              "data");

  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(), Fields(Fields::Single::kUvw)),
              "data");
  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(Fields::Single::kUvw),
                                      Fields(Fields::Single::kUvw)),
              "data, uvw");
  CheckFields(Fields(kDataField)
                  .UpdateRequirements(Fields(Fields::Single::kUvw), Fields()),
              "data, uvw");
}

BOOST_AUTO_TEST_CASE(combine) {
  Fields fields;
  fields |= Fields(Fields::Single::kFlags);
  CheckFields(fields, "flags");
  fields |= Fields(Fields::Single::kFullResFlags);
  CheckFields(fields, "flags, fullResFlags");

  const Fields data(Fields::Single::kData);
  const Fields uvw(Fields::Single::kUvw);
  CheckFields(data | uvw, "data, uvw");
}

BOOST_AUTO_TEST_SUITE_END()
