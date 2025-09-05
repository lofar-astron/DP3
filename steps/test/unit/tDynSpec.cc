// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../../common/ParameterSet.h"
#include "../../DynSpec.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::DynSpec;

BOOST_AUTO_TEST_SUITE(dynspec)

BOOST_AUTO_TEST_CASE(no_sourcelist) {
  dp3::common::ParameterSet parset;
  parset.add("sourcelist", "");

  BOOST_CHECK_THROW(std::make_unique<DynSpec>(parset, ""), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
