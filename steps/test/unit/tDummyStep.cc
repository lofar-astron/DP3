// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Template file for creating tests for a step.

#include "../../DummyStep.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::DummyStep;

BOOST_AUTO_TEST_SUITE(dummystep)

BOOST_AUTO_TEST_CASE(constructor) {
  // Test if the DummyStep can be constructed.
  // Note that this test is not a dummy: For example, if the DummyStep is
  // abstract, the compiler will complain which means the DummyStep template
  // needs updating.
  dp3::common::ParameterSet parset;
  BOOST_CHECK_NO_THROW(std::make_unique<DummyStep>(parset, "dummy"));
}

BOOST_AUTO_TEST_SUITE_END()