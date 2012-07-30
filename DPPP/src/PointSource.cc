//# PointSource.cc: Point source model component with optional spectral index
//# and rotation measure.
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
#include <DPPP/PointSource.h>
#include <DPPP/ModelComponentVisitor.h>
#include <Common/lofar_math.h>
#include <casa/BasicSL/Constants.h>

namespace LOFAR
{
namespace DPPP
{

PointSource::PointSource(const Position &position)
    :   itsPosition(position),
        itsRefFreq(0.0),
        itsPolarizedFraction(0.0),
        itsPolarizationAngle(0.0),
        itsRotationMeasure(0.0),
        itsHasRotationMeasure(false)
{
}

PointSource::PointSource(const Position &position, const Stokes &stokes)
    :   itsPosition(position),
        itsStokes(stokes),
        itsRefFreq(0.0),
        itsPolarizedFraction(0.0),
        itsPolarizationAngle(0.0),
        itsRotationMeasure(0.0),
        itsHasRotationMeasure(false)
{
}

void PointSource::setPosition(const Position &position)
{
    itsPosition = position;
}

void PointSource::setStokes(const Stokes &stokes)
{
    itsStokes = stokes;
}

void PointSource::setRotationMeasure(double fraction, double angle, double rm)
{
    itsPolarizedFraction = fraction;
    itsPolarizationAngle = angle;
    itsRotationMeasure = rm;
    itsHasRotationMeasure = true;
}

Stokes PointSource::stokes(double freq) const
{
    Stokes stokes(itsStokes);

    if(hasSpectralIndex())
    {
        // Compute spectral index as:
        // (v / v0) ^ (c0 + c1 * log10(v / v0) + c2 * log10(v / v0)^2 + ...)
        // Where v is the frequency and v0 is the reference frequency.

        // Compute log10(v / v0).
        double base = log10(freq) - log10(itsRefFreq);

        // Compute c0 + log10(v / v0) * c1 + log10(v / v0)^2 * c2 + ...
        // using Horner's rule.
        double exponent = 0.0;
        typedef vector<double>::const_reverse_iterator iterator_type;
        for(iterator_type it = itsSpectralIndex.rbegin(),
            end = itsSpectralIndex.rend(); it != end; ++it)
        {
            exponent = exponent * base + *it;
        }

        // Compute I * (v / v0) ^ exponent, where I is the value of Stokes
        // I at the reference frequency.
        stokes.I *= pow(10., base * exponent);
    }

    if(hasRotationMeasure())
    {
        double lambda = casa::C::c / freq;
        double chi = 2.0 * (itsPolarizationAngle + itsRotationMeasure
            * lambda * lambda);
        double stokesQU = stokes.I * itsPolarizedFraction;
        stokes.Q = stokesQU * cos(chi);
        stokes.U = stokesQU * sin(chi);
    }

    return stokes;
}

void PointSource::accept(ModelComponentVisitor &visitor) const
{
    visitor.visit(*this);
}

bool PointSource::hasSpectralIndex() const
{
    return !itsSpectralIndex.empty();
}

bool PointSource::hasRotationMeasure() const
{
    return itsHasRotationMeasure;
}

} //# namespace DPPP
} //# namespace LOFAR
