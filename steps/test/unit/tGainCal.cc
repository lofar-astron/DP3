// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../GainCal.h"

#include "../../../common/ParameterSet.h"
#include "mock/MockInput.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(gaincal)

BOOST_AUTO_TEST_CASE(duplicate_modelcolumn) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("gaincal.parmdb", "foo");
  parset.add("gaincal.caltype", "scalar");
  parset.add("gaincal.usemodelcolumn", "true");
  parset.add("gaincal.modelcolumn", "foo");
  BOOST_CHECK_NO_THROW(
      boost::make_unique<dp3::steps::GainCal>(input, parset, "gaincal."));

  // When the deprecated msin.modelcolumn key is also present, GainCal throws.
  parset.add("msin.modelcolumn", "bar");
  BOOST_CHECK_THROW(
      boost::make_unique<dp3::steps::GainCal>(input, parset, "gaincal."),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()