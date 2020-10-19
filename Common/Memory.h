// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Define common functions for memory (usage)
/// @author Lars Krombeen

#ifndef DP3_MEMORY_H
#define DP3_MEMORY_H

namespace DP3 {

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

}  // namespace DP3

#endif
