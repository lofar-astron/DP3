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

#include <casacore/casa/OS/HostInfo.h>

#include "Memory.h"

#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace DP3 {

double AvailableMemory(const double memory, const double memory_percentage,
                       const bool clip) {
  // Determine available memory (bytes).
  double max_system_memory = casacore::HostInfo::memoryTotal() * 1024.;

  // Determine how much memory can be used (bytes).
  double memory_max = memory * 1024 * 1024 * 1024;
  if (clip) memory_max = std::min(memory_max, max_system_memory);

  if (memory_max > max_system_memory) {
    std::cout << "WARNING: DP3 will use more memory than available."
              << std::endl;
    std::cout << max_system_memory << " bytes are available, but using "
              << memory_max << std::endl;
  }

  double memory_avail = memory_max;

  if (memory_percentage > 0) {
    memory_avail = memory_percentage * max_system_memory / 100.;
    // If the memory_avail exceeds memory_max, clip to memory_max.
    if (memory_max > 0 && memory_avail > memory_max) {
      memory_avail = memory_max;
    }
  } else if (memory <= 0) {
    // Nothing given, so use available memory on this machine.
    // Set 50% (max 2 GB) aside for other purposes.
    memory_avail = max_system_memory -
                   std::min(0.5 * max_system_memory, 2. * 1024 * 1024 * 1024);
  }

  return memory_avail;
}

}  // namespace DP3
