// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// This file contains common helper functions for BDA tests.
// tBDABuffer.cc contains the corresponding implementation, besides the
// BDABuffer tests.

#include <boost/test/unit_test.hpp>

#include "../../BDABuffer.h"

namespace dp3 {
namespace base {
namespace test {

/**
 * Verify that two BDA rows have equal metadata.
 * @param left The first BDA row in the comparison.
 * @param right The second BDA row in the comparison.
 */
void CheckBDARowMetaData(const BDABuffer::Row& left,
                         const BDABuffer::Row& right);

/**
 * Verify that all rows of two BDABuffers have equal metadata.
 * @param left The first BDABuffer in the comparison.
 * @param right The second BDABuffer in the comparison.
 */
void CheckBDARowMetaData(const BDABuffer& left, const BDABuffer& right);

}  // namespace test
}  // namespace base
}  // namespace dp3