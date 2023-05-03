// GaussianSource.h: Gaussian source model component.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_GAUSSIANSOURCE_H
#define DPPP_GAUSSIANSOURCE_H

#include "PointSource.h"

namespace dp3 {
namespace base {

/// \brief Gaussian source model component.

/// @{

class GaussianSource : public PointSource {
 public:
  typedef std::shared_ptr<GaussianSource> Ptr;
  typedef std::shared_ptr<const GaussianSource> ConstPtr;

  GaussianSource(const Direction &direction);
  GaussianSource(const Direction &direction, const Stokes &stokes);

  /// Set position angle in radians. The position angle is the smallest angle
  /// between the major axis and North, measured positively North over East.
  void setPositionAngle(double angle);
  double getPositionAngle() const { return itsPositionAngle; }

  /// Set whether the position angle (orientation) is absolute, see
  /// documentation of class member)
  void setPositionAngleIsAbsolute(bool positionAngleIsAbsolute) {
    itsPositionAngleIsAbsolute = positionAngleIsAbsolute;
  }

  /// Return whether the position angle (orientation) is absolute, see
  /// documentation of class member.
  bool getPositionAngleIsAbsolute() const { return itsPositionAngleIsAbsolute; }

  /// Set the major axis length (FWHM in radians).
  void setMajorAxis(double fwhm);
  double getMajorAxis() const { return itsMajorAxis; }

  /// Set the minor axis length (FWHM in radians).
  void setMinorAxis(double fwhm);
  double getMinorAxis() const { return itsMinorAxis; }

  void accept(ModelComponentVisitor &visitor) const override;

 private:
  double itsPositionAngle;
  /// Whether the position angle (also refered to as orientation) is absolute
  /// (w.r.t. to the local declination axis) or with respect to the declination
  /// axis at the phase center (the default until 2022, it was fixed in 5.3.0)
  bool itsPositionAngleIsAbsolute;
  double itsMajorAxis;
  double itsMinorAxis;
};

/// @}

}  // namespace base
}  // namespace dp3

#endif
