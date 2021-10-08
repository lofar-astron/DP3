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

ostream &NSTimer::print(ostream &str) const {
  if (name_.empty()) {
    str << "timer: ";
  } else {
    str << left << setw(25) << name_ << ": " << right;
  }
  if (count_ == 0) {
    str << "not used";
  } else {
    const double total = getElapsed();
    // clang-format off
    str << "avg = " << PrettyTime(total / count_)
        << ", total = " << PrettyTime(total)
        << ", count = " << setw(9) << count_;
    // clang-format on
  }
  return str;
}

}  // namespace common
}  // namespace dp3
