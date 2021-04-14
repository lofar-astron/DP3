// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// This file contains common helper functions for BDA tests.
// tBDABuffer.cc contains the corresponding implementation, besides the
// BDABuffer tests.

#include <boost/test/unit_test.hpp>

namespace dp3 {
namespace base {

class BDABuffer;

namespace test {

/**
 * Verify that the rows of two BDABuffers have equal metadata.
 * @param left The first BDABuffer in the comparison.
 * @param right The second BDABuffer in the comparison.
 */
void CheckBDARowMetaData(const BDABuffer& left, const BDABuffer& right);

}  // namespace test
}  // namespace base
}  // namespace dp3