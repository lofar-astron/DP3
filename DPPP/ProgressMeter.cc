// ProgressMeter.cc: Visual indication of a task's progress.
// Copyright (C) 1997,2000,2001,2002
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
//
// $Id$

#include "ProgressMeter.h"

#include <iostream>
#include <string>
#include <vector>

namespace DP3 {

// First implement a simple stderr based progress meter that just prints out
// 0%....10....20....30....40....50....60....70....80....90....100%
// Original design decision:
// cerr is better than cout because it isn't buffered usually, so the above
// will come out right away. Also, one often wants to direct "real" output
// to a file, but see informative messages on the screen.
// However, outputting to cerr gives quite anoying behaviour, so this
// was changed to cout in the end.

// If we have lots and lots of progress meters we should figure out
// a way to reclaim the following storage.
static std::vector<double> stderr_min, stderr_max, stderr_last;
static int stderr_creation_function(double min, double max,
    const std::string&, const std::string&,
    const std::string&, const std::string&,
    bool)
{
    stderr_min.push_back (min);
    stderr_max.push_back (max);
    stderr_last.push_back (min);
    std::cout << "\n0%" << std::flush;
    return stderr_min.size();
}

static void stderr_update_function(int id, double value)
{
    if (id < 0 || id > int(stderr_min.size())) {
      std::cerr << __FILE__ << " illegal id " << id << std::endl;
      return;
    }
    id--; // 0-relative
    int percent     = int((value - stderr_min[id]) / 
        (stderr_max[id] - stderr_min[id]) * 100.0);
    int lastpercent = int((stderr_last[id] - stderr_min[id]) / 
        (stderr_max[id] - stderr_min[id]) * 100.0);
    if (percent > lastpercent) {
      stderr_last[id] = value;
      // Probably we could do this more efficiently. We need to get all the
      // "missing" ..'s etc if we have jumped a lot since our last updated.
      for (int i=lastpercent+1; i<=percent; i++) {
        if (i%2 == 0 && i%10 != 0) {
          std::cout << "." << std::flush;
        } else if (i %10 == 0) {
          std::cout << i << std::flush;
          if (i >= 100) {
            std::cout << "%\n" << std::flush;
          }
        }
      }
    }
    
}

int (*ProgressMeter::creation_function_p)(double, double, 
    const std::string&, const std::string&,
                              const std::string&, const std::string&,
                              bool) = stderr_creation_function;
void (*ProgressMeter::update_function_p)(int, double) = stderr_update_function;

ProgressMeter::ProgressMeter()
    : id_p(-1), min_p(0.0), max_p(1.0), update_every_p(1), update_count_p(0)
{
}

ProgressMeter::ProgressMeter(double min, double max, 
    const std::string& title, const std::string& subtitle,
    const std::string& minlabel, const std::string& maxlabel,
    bool estimateTime, int updateEvery)
    : id_p(-1), min_p(min), max_p(max), update_every_p(updateEvery),
      update_count_p(0)
{
    // Correct silently
    if (update_every_p <= 0) {
      update_every_p = 1;
    }
    if (creation_function_p) {
      id_p = creation_function_p(min, max, title, subtitle, minlabel, maxlabel,
          estimateTime);
    }
}

ProgressMeter::~ProgressMeter()
{
    update_count_p++;
    update(max_p, true);
}

void ProgressMeter::update(double value, bool force)
{
    update_count_p++;
    // Always force the first one through
    if (update_count_p == 1) {
      force = true;
    }
    if((value >= min_p) && (value <= max_p)){
      if(update_count_p == 1 || force || ((update_count_p%update_every_p)== 0))
      {
        // Do the update if we have a "sink" and a valid id
        if (id_p > 0 && update_function_p) {
          update_function_p(id_p, value);
        } else {
          // If we have more than one progress meter active at once
          // this might look pretty confusing. We can decide what to
          // do if that ever actually happens.

        }
      }
    }
    else{

      //cerr << "WARNING: progress meter trying to update beyond range" << endl;//The user does not need to know that the programmer does not know how to add.
      
    }
}

} // NAMESPACE LOFAR - END
