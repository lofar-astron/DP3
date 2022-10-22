// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Split.h"

#include <boost/test/unit_test.hpp>

#include <dp3/base/DPInfo.h>

#include "../../../common/ParameterSet.h"

#include "../../NullStep.h"
#include "../../PhaseShift.h"
#include "../../Upsample.h"
#include "mock/MockInput.h"

using dp3::steps::Split;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(split)

// Test that split works when the sub-step list is empty.
BOOST_AUTO_TEST_CASE(no_sub_steps) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("split.replaceparms", "foo_param");
  parset.add("foo_param", "[foo,bar]");
  parset.add("split.steps", "[]");
  Split no_steps(&input, parset, "split.");

  BOOST_CHECK_NO_THROW(no_steps.updateInfo(dp3::base::DPInfo()));
  BOOST_CHECK_NO_THROW(no_steps.process(dp3::base::DPBuffer()));
  BOOST_CHECK_NO_THROW(no_steps.finish());
}

BOOST_AUTO_TEST_CASE(fields) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;

  // Create some steps with different non-empty required fields.
  parset.add("phaseshift.phasecenter", "foo_center");
  parset.add("upsample.timestep", "2");
  parset.add("upsample.updateuvw", "true");
  const dp3::steps::PhaseShift phase_shift(parset, "phaseshift.");
  const dp3::steps::Upsample upsample(parset, "upsample.");
  BOOST_TEST(phase_shift.getRequiredFields() != dp3::common::Fields());
  BOOST_TEST(upsample.getRequiredFields() != dp3::common::Fields());
  BOOST_TEST(phase_shift.getRequiredFields() != upsample.getRequiredFields());

  parset.add("split.replaceparms", "foo_param");
  parset.add("foo_param", "[foo,bar]");
  parset.add("split.steps", "[]");
  const Split no_steps(&input, parset, "split.");
  BOOST_TEST(no_steps.getRequiredFields() == dp3::common::Fields());

  parset.replace("split.steps", "[phaseshift]");
  const Split one_step(&input, parset, "split.");
  BOOST_TEST(one_step.getRequiredFields() == phase_shift.getRequiredFields());

  parset.replace("split.steps", "[phaseshift,upsample]");
  const Split two_steps(&input, parset, "split.");
  BOOST_TEST(two_steps.getRequiredFields() ==
             (phase_shift.getRequiredFields() | upsample.getRequiredFields()));
}

BOOST_AUTO_TEST_SUITE_END()
