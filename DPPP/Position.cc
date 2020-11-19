// Position.cc: A position on the celestial sphere.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "Position.h"

namespace DP3 {
namespace DPPP {

Position::Position() {
  itsPosition[0] = 0.0;
  itsPosition[1] = 0.0;
}

Position::Position(double alpha, double delta) {
  itsPosition[0] = alpha;
  itsPosition[1] = delta;
}

}  // namespace DPPP
}  // namespace DP3
