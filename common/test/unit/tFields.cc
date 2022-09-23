// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Fields.h"

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

BOOST_AUTO_TEST_SUITE_END()
