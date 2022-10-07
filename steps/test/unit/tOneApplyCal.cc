// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../OneApplyCal.h"

#include <boost/test/unit_test.hpp>

#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

using dp3::steps::OneApplyCal;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(one_apply_cal)

BOOST_AUTO_TEST_CASE(fields) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("parmdb", "tApplyCal_tmp.parmdb");

  const OneApplyCal one_apply_cal(&input, parset, "", "");
  BOOST_TEST(one_apply_cal.getRequiredFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
  BOOST_TEST(one_apply_cal.getProvidedFields() ==
             (Step::kDataField | Step::kFlagsField));

  parset.add("updateweights", "true");
  const OneApplyCal updates_weights(&input, parset, "", "");
  BOOST_TEST(updates_weights.getRequiredFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
  BOOST_TEST(updates_weights.getProvidedFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
}

BOOST_AUTO_TEST_SUITE_END()
