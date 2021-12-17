// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_PIERCEPOINT_H_
#define DP3_DDECAL_PIERCEPOINT_H_

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <armadillo>

#include <vector>

namespace dp3 {
namespace ddecal {

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
  double getHeight() const { return itsIonoHeight; }

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

}  // namespace ddecal
}  // namespace dp3

#endif
