// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This file has generic helper code for tests.

#ifndef DP3_COMMON_TEST_UNIT_TCOMMON_H_
#define DP3_COMMON_TEST_UNIT_TCOMMON_H_

#include <boost/test/data/test_case.hpp>

namespace dp3::common::test {

/// Data range used in many tests.
inline const auto kTrueFalseRange = boost::unit_test::data::make({true, false});

}  // namespace dp3::common::test

#endif
