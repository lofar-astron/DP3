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

#ifndef PIERCEPOINT_H
#define PIERCEPOINT_H

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <armadillo>

#include <vector>

namespace DP3 {

class PiercePoint {
  /// default height in meter (300000)
  static const double IONOheight;

  /// default Earth radius in meter (6371000)
  static const double EarthRadius;

 public:
  PiercePoint(double height = PiercePoint::IONOheight);
  PiercePoint(const casacore::MPosition &ant,
              const casacore::MDirection &source, const double height);
  PiercePoint(const casacore::MPosition &ant,
              const casacore::MDirection &source);

  void init(const casacore::MPosition &ant, const casacore::MDirection &source,
            const double height);

  void evaluate(casacore::MEpoch time);
  arma::Col<double> getValue() const { return itsValue; }
  casacore::MPosition getPos() const { return itsPosition; }
  casacore::MDirection getDir() const { return itsDirection; }

 private:
  /// Station position.
  casacore::MPosition itsPosition;
  /// Source position.
  casacore::MDirection itsDirection;
  /// Ionospheric layer height.
  double itsIonoHeight;
  /// Square of length antenna vector (int ITRF) minus square of vector to
  /// piercepoint. This is constant for an assumed spherical Earth.
  double itsC;
  arma::Col<double> itsValue;  ///< PiercePoint in ITRF coordinates
};

}  // namespace DP3

#endif
