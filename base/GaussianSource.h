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
  double positionAngle() const;

  /// Set the major axis length (FWHM in radians).
  void setMajorAxis(double fwhm);
  double majorAxis() const;

  /// Set the minor axis length (FWHM in radians).
  void setMinorAxis(double fwhm);
  double minorAxis() const;

  virtual void accept(ModelComponentVisitor &visitor) const;

 private:
  double itsPositionAngle;
  double itsMajorAxis;
  double itsMinorAxis;
};

/// @}

// -------------------------------------------------------------------------- //
// - Implementation: GaussianSource                                         - //
// -------------------------------------------------------------------------- //

inline double GaussianSource::positionAngle() const { return itsPositionAngle; }

inline double GaussianSource::majorAxis() const { return itsMajorAxis; }

inline double GaussianSource::minorAxis() const { return itsMinorAxis; }

}  // namespace base
}  // namespace dp3

#endif
