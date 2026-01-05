// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Patch.h"

#include "base/ModelComponentVisitor.h"

#include <cmath>

namespace dp3::model {

void Patch::ComputeDirection() {
  if (components_.empty()) {
    direction_ = base::Direction();
  } else {
    double x = 0.0, y = 0.0, z = 0.0;
    for (const auto& component : components_) {
      const base::Direction& direction = component->direction();
      const double cos_dec = std::cos(direction.dec);
      x += std::cos(direction.ra) * cos_dec;
      y += std::sin(direction.ra) * cos_dec;
      z += std::sin(direction.dec);
    }

    x /= components_.size();
    y /= components_.size();
    z /= components_.size();

    direction_ = base::Direction(std::atan2(y, x), std::asin(z));
  }
}

}  // namespace dp3::model
