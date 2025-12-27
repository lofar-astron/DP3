// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// This file contains common helper functions for BDA tests.
// tBdaBuffer.cc contains the corresponding implementation, besides the
// BdaBuffer tests.

#include <boost/test/unit_test.hpp>

#include "base/BdaBuffer.h"

namespace dp3 {
namespace base {
namespace test {

/**
 * Verify that two BDA rows have equal metadata.
 * @param left The first BDA row in the comparison.
 * @param right The second BDA row in the comparison.
 */
void CheckBDARowMetaData(const BdaBuffer::Row& left,
                         const BdaBuffer::Row& right);

/**
 * Verify that all rows of two BdaBuffers have equal metadata.
 * @param left The first BdaBuffer in the comparison.
 * @param right The second BdaBuffer in the comparison.
 */
void CheckBDARowMetaData(const BdaBuffer& left, const BdaBuffer& right);

}  // namespace test
}  // namespace base
}  // namespace dp3
