// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../ColumnReader.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

using dp3::steps::ColumnReader;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(column_reader)

BOOST_AUTO_TEST_CASE(fields_default) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  const ColumnReader column_reader(input, parset, "");
  BOOST_TEST(column_reader.getRequiredFields() == dp3::common::Fields());
  BOOST_TEST(column_reader.getProvidedFields() == Step::kDataField);
}

BOOST_DATA_TEST_CASE(fields_add_subtract,
                     boost::unit_test::data::make({"add", "subtract"}),
                     operation) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("operation", operation);
  const ColumnReader column_reader(input, parset, "");
  BOOST_TEST(column_reader.getRequiredFields() == Step::kDataField);
  BOOST_TEST(column_reader.getProvidedFields() == Step::kDataField);
}

BOOST_AUTO_TEST_SUITE_END()
