// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Memory.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::common::AvailableMemory;

namespace {
constexpr double kGB2BFactor = 1024 * 1024 * 1024;
constexpr double kTooMuchMemory = 1000 * 1024;                    // GB
constexpr double kTooMuchMemoryB = kTooMuchMemory * kGB2BFactor;  // B
}  // namespace

BOOST_AUTO_TEST_SUITE(memory)

// By default the returned value should be clipped and therefore be smaller than
// the input.
BOOST_AUTO_TEST_CASE(clipping) {
  double mem = AvailableMemory(kTooMuchMemory);
  BOOST_TEST(mem < kTooMuchMemoryB);
}

// When not clipping, the result can be equal to the input.
BOOST_AUTO_TEST_CASE(no_clipping) {
  double mem = AvailableMemory(kTooMuchMemory, 0, false);
  BOOST_TEST(mem == kTooMuchMemoryB);
}

// When the percentage memory result exceeds gbs the returned result should be
// gbs.
BOOST_AUTO_TEST_CASE(percentage) {
  double bytes = 200;
  double gbs = bytes / kGB2BFactor;
  double mem = AvailableMemory(gbs, 80);
  BOOST_TEST(mem == bytes);
}

// When no parameters are passed, the function should return approximately 50%
// of the available memory.
BOOST_AUTO_TEST_CASE(default_test) {
  double max_mem = AvailableMemory(kTooMuchMemory);
  double default_mem = AvailableMemory();

  double expected_mem = max_mem - std::min(0.5 * max_mem, 2. * kGB2BFactor);

  BOOST_TEST(default_mem == expected_mem);
}

BOOST_AUTO_TEST_CASE(mem_too_small) {
  BOOST_CHECK_THROW(AvailableMemory(0, -0.01), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(mem_too_large) {
  BOOST_CHECK_THROW(AvailableMemory(0, 100.01), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
