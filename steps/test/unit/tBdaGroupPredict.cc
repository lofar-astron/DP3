// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../BdaGroupPredict.h"

#include <boost/test/unit_test.hpp>

#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

using dp3::steps::BdaGroupPredict;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(bda_group_predict)

BOOST_AUTO_TEST_CASE(fields) {
  dp3::common::ParameterSet parset;
  dp3::steps::MockInput input;
  const BdaGroupPredict bda_group_predict(input, parset, "");
  BOOST_TEST(bda_group_predict.getRequiredFields() == dp3::common::Fields());
  BOOST_TEST(bda_group_predict.getProvidedFields() == Step::kDataField);
}

BOOST_AUTO_TEST_SUITE_END()
