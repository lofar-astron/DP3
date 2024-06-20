// ProgressMeter.cc: Visual indication of a task's progress.
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "ProgressMeter.h"

#include <iostream>
#include <string>
#include <vector>

#include <aocommon/logger.h>

using aocommon::Logger;

namespace dp3 {
namespace base {

// If we have lots and lots of progress meters we should figure out
// a way to reclaim the following storage.
static std::vector<double> stderr_min, stderr_max, stderr_last;
static int stderr_creation_function(double min, double max, const std::string&,
                                    const std::string&, const std::string&,
                                    const std::string&, bool) {
  stderr_min.push_back(min);
  stderr_max.push_back(max);
  stderr_last.push_back(min);
  Logger::Info << "\n0%";
  Logger::Info.Flush();
  return stderr_min.size();
}

static void stderr_update_function(int id, double value) {
  if (id < 0 || id > int(stderr_min.size())) {
    Logger::Error << __FILE__ << " illegal id " << id << '\n';
    return;
  }
  id--;  // 0-relative
  int percent =
      int((value - stderr_min[id]) / (stderr_max[id] - stderr_min[id]) * 100.0);
  int lastpercent = int((stderr_last[id] - stderr_min[id]) /
                        (stderr_max[id] - stderr_min[id]) * 100.0);
  if (percent > lastpercent) {
    stderr_last[id] = value;
    // Probably we could do this more efficiently. We need to get all the
    // "missing" ..'s etc if we have jumped a lot since our last updated.
    for (int i = lastpercent + 1; i <= percent; i++) {
      if (i % 2 == 0 && i % 10 != 0) {
        Logger::Info << ".";
        Logger::Info.Flush();
      } else if (i % 10 == 0) {
        Logger::Info << i;
        Logger::Info.Flush();
        if (i >= 100) {
          Logger::Info << "%\n";
          Logger::Info.Flush();
        }
      }
    }
  }
}

int (*ProgressMeter::creation_function_p)(double, double, const std::string&,
                                          const std::string&,
                                          const std::string&,
                                          const std::string&,
                                          bool) = stderr_creation_function;
void (*ProgressMeter::update_function_p)(int, double) = stderr_update_function;

ProgressMeter::ProgressMeter()
    : id_p(-1), min_p(0.0), max_p(1.0), update_every_p(1), update_count_p(0) {}

ProgressMeter::ProgressMeter(double min, double max, const std::string& title,
                             const std::string& subtitle,
                             const std::string& minlabel,
                             const std::string& maxlabel, bool estimateTime,
                             int updateEvery)
    : id_p(-1),
      min_p(min),
      max_p(max),
      update_every_p(updateEvery),
      update_count_p(0) {
  // Correct silently
  if (update_every_p <= 0) {
    update_every_p = 1;
  }
  if (creation_function_p) {
    id_p = creation_function_p(min, max, title, subtitle, minlabel, maxlabel,
                               estimateTime);
  }
}

ProgressMeter::~ProgressMeter() {
  update_count_p++;
  update(max_p, true);
}

void ProgressMeter::update(double value, bool force) {
  update_count_p++;
  // Always force the first one through
  if (update_count_p == 1) {
    force = true;
  }
  if ((value >= min_p) && (value <= max_p)) {
    if (update_count_p == 1 || force ||
        ((update_count_p % update_every_p) == 0)) {
      // Do the update if we have a "sink" and a valid id
      if (id_p > 0 && update_function_p) {
        update_function_p(id_p, value);
      } else {
        // If we have more than one progress meter active at once
        // this might look pretty confusing. We can decide what to
        // do if that ever actually happens.
      }
    }
  } else {
  }
}

}  // namespace base
}  // namespace dp3
