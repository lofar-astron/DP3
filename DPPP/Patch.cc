// Patch.cc: A set of sources for which direction dependent effects are assumed
// to be equal.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "Patch.h"
#include "ModelComponentVisitor.h"

#include <cmath>

namespace DP3 {
namespace DPPP {

void Patch::computePosition() {
  itsPosition = Position();

  if (!itsComponents.empty()) {
    double x = 0.0, y = 0.0, z = 0.0;
    for (const_iterator it = begin(), it_end = end(); it != it_end; ++it) {
      const Position &position = (*it)->position();
      double cosDec = cos(position[1]);
      x += cos(position[0]) * cosDec;
      y += sin(position[0]) * cosDec;
      z += sin(position[1]);
    }

    x /= itsComponents.size();
    y /= itsComponents.size();
    z /= itsComponents.size();

    itsPosition[0] = atan2(y, x);
    itsPosition[1] = asin(z);
  }
}

}  // namespace DPPP
}  // namespace DP3
