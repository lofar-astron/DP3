// Copyright (C) 2020
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

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"

#include "DPInfo.h"
#include "SetBeam.h"

#include <casacore/casa/Quanta/MVAngle.h>

namespace DP3 {
namespace DPPP {

SetBeam::SetBeam(DPInput* input, const ParameterSet& parset,
                 const string& prefix)
    : _input(input),
      _name(prefix),
      _directionStr(parset.getStringVector(prefix + "direction",
                                           std::vector<std::string>())),
      _mode(StringToBeamCorrectionMode(
          parset.getString(prefix + "beammode", "default"))) {}

void SetBeam::updateInfo(const DPInfo& dpInfo) {
  info() = dpInfo;

  // Parse direction parset value
  if (_directionStr.empty())
    _direction = info().phaseCenter();
  else {
    if (_directionStr.size() != 2)
      throw std::runtime_error(
          "2 values must be given in direction option of SetBeam");
    casacore::MDirection phaseCenter;
    casacore::Quantity q0, q1;
    if (!casacore::MVAngle::read(q0, _directionStr[0]))
      throw std::runtime_error(
          _directionStr[0] +
          " is an invalid RA or longitude in SetBeam direction");
    if (!casacore::MVAngle::read(q1, _directionStr[1]))
      throw std::runtime_error(
          _directionStr[1] +
          " is an invalid DEC or latitude in SetBeam direction");
    casacore::MDirection::Types type = casacore::MDirection::J2000;
    _direction = casacore::MDirection(q0, q1, type);
  }

  info().setBeamCorrectionMode(_mode);
  info().setBeamCorrectionDir(_direction);
}

void SetBeam::show(std::ostream& os) const {
  os << "SetBeam " << _name << '\n'
     << "  mode:              " << BeamCorrectionModeToString(_mode) << '\n'
     << "  direction:         " << _directionStr << '\n';
}

bool SetBeam::process(const DPBuffer& buffer) {
  getNextStep()->process(buffer);
  return false;
}

}  // namespace DPPP
}  // namespace DP3
