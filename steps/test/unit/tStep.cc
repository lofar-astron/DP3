// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/Step.h"

#include <boost/test/unit_test.hpp>

#include "steps/test/unit/mock/ThrowStep.h"

using dp3::steps::test::ThrowStep;

BOOST_AUTO_TEST_SUITE(step)

BOOST_AUTO_TEST_CASE(prev_next_pointers) {
  auto step1 = std::make_shared<ThrowStep>();
  auto step2 = std::make_shared<ThrowStep>();
  BOOST_CHECK(!step1->getNextStep());
  BOOST_CHECK(!step1->getPrevStep());
  BOOST_CHECK(!step2->getNextStep());
  BOOST_CHECK(!step2->getPrevStep());

  step1->setNextStep(step2);
  BOOST_CHECK_EQUAL(step1->getNextStep(), step2);
  BOOST_CHECK(!step1->getPrevStep());
  BOOST_CHECK(!step2->getNextStep());
  BOOST_CHECK_EQUAL(step2->getPrevStep(), step1.get());

  // The destructor for step1 should update step2's previous step pointer.
  step1.reset();
  BOOST_CHECK(!step2->getPrevStep());
  BOOST_CHECK(!step2->getNextStep());
}

BOOST_AUTO_TEST_SUITE_END()
