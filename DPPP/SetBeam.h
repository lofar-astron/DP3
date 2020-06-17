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

/// @file
/// @brief DPPP step class to set the beam keywords in a ms

#ifndef DPPP_SETBEAM_H
#define DPPP_SETBEAM_H

#include "DPInput.h"
#include "DPBuffer.h"
#include "Position.h"

#include <casacore/measures/Measures/MDirection.h>

namespace DP3 {

class ParameterSet;

namespace DPPP {

/// @brief DPPP step class to set the beam keywords in a ms
class SetBeam final : public DPStep
{
public:
  /// Parameters are obtained from the parset using the given prefix.
  SetBeam (DPInput* input, const ParameterSet& parameters, const string& prefix);

  bool process(const DPBuffer& buffer) override;

  void finish() override { };

  void updateInfo(const DPInfo& info) override;

  void show(std::ostream&) const override;
  
private:
  DPInput* _input;
  string _name;
  std::vector<string> _directionStr;
  casacore::MDirection _direction;
  BeamCorrectionMode _mode;
};

} } // end namespaces

#endif

