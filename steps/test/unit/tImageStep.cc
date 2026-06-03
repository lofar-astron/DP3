// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../../common/ParameterSet.h"
#include "../../ImageStep.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::ImageStep;

namespace {
const float kCadence{4};
const std::string kOptions{"-size 1024 1024 -scale 1amin"};
}  // namespace

BOOST_AUTO_TEST_SUITE(imagestep)

BOOST_AUTO_TEST_CASE(constructor) {
  dp3::common::ParameterSet parset;

  parset.add("cadence", std::to_string(kCadence));
  parset.add("options", kOptions);

  BOOST_CHECK_NO_THROW(std::make_unique<ImageStep>(parset, ""));
}

BOOST_AUTO_TEST_SUITE_END()
