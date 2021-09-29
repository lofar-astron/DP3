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

namespace dp3 {
namespace base {

void Patch::computeDirection() {
  itsDirection = Direction();

  if (!itsComponents.empty()) {
    double x = 0.0, y = 0.0, z = 0.0;
    for (const auto& component : itsComponents) {
      const Direction& direction = component->direction();
      const double cos_dec = cos(direction.dec);
      x += cos(direction.ra) * cos_dec;
      y += sin(direction.ra) * cos_dec;
      z += sin(direction.dec);
    }

    x /= itsComponents.size();
    y /= itsComponents.size();
    z /= itsComponents.size();

    itsDirection.ra = atan2(y, x);
    itsDirection.dec = asin(z);
  }
}

}  // namespace base
}  // namespace dp3
