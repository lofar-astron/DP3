//# GaussianSource.cc: Gaussian source model component.
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

#include <lofar_config.h>
#include <DPPP/GaussianSource.h>
#include <DPPP/ModelComponentVisitor.h>

namespace LOFAR
{
namespace DPPP
{

GaussianSource::GaussianSource(const Position &position)
    :   PointSource(position),
        itsPositionAngle(0.0),
        itsMajorAxis(0.0),
        itsMinorAxis(0.0)
{
}

GaussianSource::GaussianSource(const Position &position, const Stokes &stokes)
    :   PointSource(position, stokes),
        itsPositionAngle(0.0),
        itsMajorAxis(0.0),
        itsMinorAxis(0.0)
{
}

void GaussianSource::setPositionAngle(double angle)
{
    itsPositionAngle = angle;
}

void GaussianSource::setMajorAxis(double fwhm)
{
    itsMajorAxis = fwhm;
}

void GaussianSource::setMinorAxis(double fwhm)
{
    itsMinorAxis = fwhm;
}

void GaussianSource::accept(ModelComponentVisitor &visitor) const
{
    visitor.visit(*this);
}

} //# namespace DPPP
} //# namespace LOFAR
