// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SetBeam.h"

#include <iostream>

#include <casacore/casa/Quanta/MVAngle.h>

#include "base/DPInfo.h"

#include "common/ParameterSet.h"
#include "common/StreamUtil.h"

using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

SetBeam::SetBeam(const common::ParameterSet& parset, const std::string& prefix)
    : name_(prefix),
      direction_strings_(parset.getStringVector(prefix + "direction",
                                                std::vector<std::string>())),
      mode_(everybeam::ParseBeamMode(
          parset.getString(prefix + "beammode", "default"))) {}

void SetBeam::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);

  // Parse direction parset value
  if (direction_strings_.empty())
    direction_ = getInfoOut().phaseCenter();
  else {
    if (direction_strings_.size() != 2)
      throw std::runtime_error(
          "2 values must be given in direction option of SetBeam");
    casacore::MDirection phaseCenter;
    casacore::Quantity q0, q1;
    if (!casacore::MVAngle::read(q0, direction_strings_[0]))
      throw std::runtime_error(
          direction_strings_[0] +
          " is an invalid RA or longitude in SetBeam direction");
    if (!casacore::MVAngle::read(q1, direction_strings_[1]))
      throw std::runtime_error(
          direction_strings_[1] +
          " is an invalid DEC or latitude in SetBeam direction");
    casacore::MDirection::Types type = casacore::MDirection::J2000;
    direction_ = casacore::MDirection(q0, q1, type);
  }

  GetWritableInfoOut().setBeamCorrectionMode(static_cast<int>(mode_));
  GetWritableInfoOut().setBeamCorrectionDir(direction_);
}

void SetBeam::show(std::ostream& os) const {
  os << "SetBeam " << name_ << '\n'
     << "  mode:              " << everybeam::ToString(mode_) << '\n'
     << "  direction:         " << direction_strings_ << '\n';
}

bool SetBeam::process(std::unique_ptr<dp3::base::DPBuffer> buffer) {
  getNextStep()->process(std::move(buffer));
  return false;
}

}  // namespace steps
}  // namespace dp3
