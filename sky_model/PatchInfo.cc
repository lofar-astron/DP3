// PatchInfo.cc: Info about a patch
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// @file
// @brief Info about a patch
// @author Ger van Diepen (diepen AT astron nl)

#include "PatchInfo.h"

#include <casacore/casa/Quanta/MVAngle.h>

#include <cassert>

using casacore::MVAngle;

namespace dp3::sky_model {

void PatchSumInfo::add(double ra, double dec, double flux) {
  // Add the position in xyz coordinates.
  // Use the flux as the weight.
  double cosDec = std::cos(dec);
  itsSumX += std::cos(ra) * cosDec * flux;
  itsSumY += std::sin(ra) * cosDec * flux;
  itsSumZ += std::sin(dec) * flux;
  itsSumFlux += flux;
}

void toSkyModel(std::ostream& output, const PatchInfo& patch) {
  output << ", , " << patch.getName() << ", ";
  MVAngle(patch.getRa()).print(output, MVAngle::Format(MVAngle::TIME, 9));
  output << ", ";
  MVAngle(patch.getDec()).print(output, MVAngle::Format(MVAngle::ANGLE, 9));
  output << '\n';
}

}  // namespace dp3::sky_model
