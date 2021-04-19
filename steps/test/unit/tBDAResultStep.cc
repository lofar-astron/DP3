// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../BDAResultStep.h"
#include "../../../base/BDABuffer.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

using dp3::base::BDABuffer;
using dp3::steps::BDAResultStep;

BOOST_AUTO_TEST_SUITE(bdaresultstep)

BOOST_AUTO_TEST_CASE(constructor) {
  BDAResultStep step;
  BOOST_CHECK(step.Extract().empty());
}

BOOST_AUTO_TEST_CASE(process_and_extract) {
  BDAResultStep step;

  auto buffer1 = boost::make_unique<BDABuffer>(0, BDABuffer::Fields(false));
  auto buffer2 = boost::make_unique<BDABuffer>(0, BDABuffer::Fields(false));
  BDABuffer* const buffer1_ptr = buffer1.get();
  BDABuffer* const buffer2_ptr = buffer2.get();

  BOOST_CHECK(step.process(std::move(buffer1)));
  BOOST_CHECK(!buffer1);
  BOOST_CHECK(step.process(std::move(buffer2)));
  BOOST_CHECK(!buffer2);

  // Check that Extract() returns the buffers in the right order.
  std::vector<std::unique_ptr<BDABuffer>> extracted = step.Extract();
  BOOST_REQUIRE(extracted.size() == 2U);
  BOOST_CHECK(extracted[0].get() == buffer1_ptr);
  BOOST_CHECK(extracted[1].get() == buffer2_ptr);

  BOOST_CHECK(step.Extract().empty());
}

BOOST_AUTO_TEST_SUITE_END()