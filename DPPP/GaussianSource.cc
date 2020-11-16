// GaussianSource.cc: Gaussian source model component.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "GaussianSource.h"
#include "ModelComponentVisitor.h"

namespace DP3 {
namespace DPPP {

GaussianSource::GaussianSource(const Position &position)
    : PointSource(position),
      itsPositionAngle(0.0),
      itsMajorAxis(0.0),
      itsMinorAxis(0.0) {}

GaussianSource::GaussianSource(const Position &position, const Stokes &stokes)
    : PointSource(position, stokes),
      itsPositionAngle(0.0),
      itsMajorAxis(0.0),
      itsMinorAxis(0.0) {}

void GaussianSource::setPositionAngle(double angle) {
  itsPositionAngle = angle;
}

void GaussianSource::setMajorAxis(double fwhm) { itsMajorAxis = fwhm; }

void GaussianSource::setMinorAxis(double fwhm) { itsMinorAxis = fwhm; }

void GaussianSource::accept(ModelComponentVisitor &visitor) const {
  visitor.visit(*this);
}

}  // namespace DPPP
}  // namespace DP3
