//# GaussianSource.h: Gaussian source model component.
//#
//# Copyright (C) 2012
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#ifndef DPPP_GAUSSIANSOURCE_H
#define DPPP_GAUSSIANSOURCE_H

// \file
// Gaussian source model component.

#include <DPPP/PointSource.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

class GaussianSource: public PointSource
{
public:
    typedef shared_ptr<GaussianSource>          Ptr;
    typedef shared_ptr<const GaussianSource>    ConstPtr;

    GaussianSource(const Position &position);
    GaussianSource(const Position &position, const Stokes &stokes);

    // Set position angle in radians. The position angle is the smallest angle
    // between the major axis and North, measured positively North over East.
    void setPositionAngle(double angle);
    double positionAngle() const;

    // Set the major axis length (FWHM in radians).
    void setMajorAxis(double fwhm);
    double majorAxis() const;

    // Set the minor axis length (FWHM in radians).
    void setMinorAxis(double fwhm);
    double minorAxis() const;

    virtual void accept(ModelComponentVisitor &visitor) const;

private:
    double      itsPositionAngle;
    double      itsMajorAxis;
    double      itsMinorAxis;
};

// @}

// -------------------------------------------------------------------------- //
// - Implementation: GaussianSource                                         - //
// -------------------------------------------------------------------------- //

inline double GaussianSource::positionAngle() const
{
    return itsPositionAngle;
}

inline double GaussianSource::majorAxis() const
{
    return itsMajorAxis;
}

inline double GaussianSource::minorAxis() const
{
    return itsMinorAxis;
}

} //# namespace DPPP
} //# namespace LOFAR

#endif
