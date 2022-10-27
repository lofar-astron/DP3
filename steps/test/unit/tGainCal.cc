// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../GainCal.h"

#include <boost/test/unit_test.hpp>

#include "../../../common/ParameterSet.h"

using dp3::steps::GainCal;

BOOST_AUTO_TEST_SUITE(gaincal)

BOOST_AUTO_TEST_CASE(duplicate_modelcolumn) {
  dp3::common::ParameterSet parset;
  parset.add("gaincal.parmdb", "foo");
  parset.add("gaincal.caltype", "scalar");
  parset.add("gaincal.usemodelcolumn", "true");
  parset.add("gaincal.modelcolumn", "foo");
  BOOST_CHECK_NO_THROW(std::make_unique<GainCal>(parset, "gaincal."));

  // When the deprecated msin.modelcolumn key is also present, GainCal throws.
  parset.add("msin.modelcolumn", "bar");
  BOOST_CHECK_THROW(std::make_unique<GainCal>(parset, "gaincal."),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(provided_fields) {
  dp3::common::ParameterSet parset;
  parset.add("parmdb", "foo");
  parset.add("caltype", "scalar");
  parset.add("usemodelcolumn", "true");
  parset.add("modelcolumn", "bar");
  const GainCal gaincal(parset, "");
  BOOST_TEST(gaincal.getProvidedFields() == dp3::common::Fields());

  parset.add("applysolution", "true");
  const GainCal applies_solution(parset, "");
  BOOST_TEST(applies_solution.getProvidedFields() ==
             (GainCal::kDataField | GainCal::kFlagsField));
}

BOOST_AUTO_TEST_SUITE_END()
