// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Define common functions for memory (usage)
/// @author Lars Krombeen

#ifndef DP3_MEMORY_H
#define DP3_MEMORY_H

namespace dp3 {
namespace common {

/**
 * Determines the available memory that can be used by the program.
 * It takes two (optional) inputs.
 *
 * If memory_percentage the result will be the maximum amount of
 * available memory of the system * percentage / 100. It will be clipped to
 * memory if set.
 *
 * If memory is set, the system that value will be returned, unless it is
 * greater than the available memory of the system in which case that number
 * will be returned.
 *
 * If no parameters are set, the returned value will be the available system
 * memoery minus a small margin of max 2GB.
 *
 * @param memory, amount of wanted memory in GB.
 * @param memory_percentage, >= 0 && <= 100, percentage of the available system
 * memory to use.
 * @param clip, if true, no more memory will be used than available on the
 * system. Limit can still be exceeded if multiple steps use 100% for example.
 * @return Recommended amount of memory to use in BYTES.
 */
double AvailableMemory(const double memory = 0,
                       const double memory_percentage = 0,
                       const bool clip = true);

}  // namespace common
}  // namespace dp3

#endif
