// Timer.cc: Accurate timer
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "PrettyUnits.h"
#include "Timer.h"

using namespace std;

namespace dp3 {
namespace common {

double NSTimer::CPU_speed_in_MHz = NSTimer::get_CPU_speed_in_MHz();

double NSTimer::get_CPU_speed_in_MHz() {
  // first a few sanity checks
  static_assert(sizeof(int) == 4, "sizeof(int) == 4 is required");
  static_assert(sizeof(long long) == 8, "sizeof(long long) == 8 is required");

#if (defined __linux__ || defined __blrts__) &&                               \
    (defined __i386__ || defined __x86_64__ || defined __ia64__ ||            \
     defined __PPC__) &&                                                      \
    (defined __GNUC__ || defined __INTEL_COMPILER || defined __PATHSCALE__ || \
     defined __xlC__)
  ifstream infile("/proc/cpuinfo");
  char buffer[256], *colon;

  while (infile.good()) {
    infile.getline(buffer, sizeof buffer);

#if defined __PPC__
    if (strcmp("cpu\t\t: 450 Blue Gene/P DD2", buffer) == 0) return 850.0;

    if (strcmp("machine\t\t: Blue Gene", buffer) == 0) return 700.0;

    if (strncmp("timebase", buffer, 8) == 0 &&
        (colon = strchr(buffer, ':')) != 0)
      return atof(colon + 2) / 1e6;
#else
    if (strncmp("cpu MHz", buffer, 7) == 0 &&
        (colon = strchr(buffer, ':')) != 0)
      return atof(colon + 2);
#endif
  }

  return 0.0;
#else
#warning unsupported architecture
  return 0.0;
#endif
}

double NSTimer::getElapsed() const {
  double time = total_time.full / 1e6;
  if (CPU_speed_in_MHz > 0) {
    time /= CPU_speed_in_MHz;
  }
  return time;
}

ostream &NSTimer::print(ostream &str) const {
  if (itsName.empty()) {
    str << "timer: ";
  } else {
    str << left << setw(25) << itsName << ": " << right;
  }
  if (count == 0) {
    str << "not used";
  } else {
    double total = static_cast<double>(total_time.full);
    if (CPU_speed_in_MHz == 0) {
      str << "avg = " << total / static_cast<double>(count);
      str << ", total = " << total_time.full << " cycles";
    } else {
      total /= 1e6 * CPU_speed_in_MHz;
      str << "avg = " << PrettyTime(total / static_cast<double>(count))
          << ", total = " << PrettyTime(total);
    }
    str << ", count = " << setw(9) << count;
  }
  return str;
}

}  // namespace common
}  // namespace dp3
