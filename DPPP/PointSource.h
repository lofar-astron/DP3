//# PointSource.h: Point source model component with optional spectral index and
//# rotation measure.
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

#ifndef DPPP_POINTSOURCE_H
#define DPPP_POINTSOURCE_H

#include <vector>

// \file
// Point source model component with optional spectral index and rotation
// measure.

#include "ModelComponent.h"
#include "Position.h"
#include "Stokes.h"

#include <memory>

namespace DP3
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

class PointSource: public ModelComponent
{
public:
    typedef std::shared_ptr<PointSource>         Ptr;
    typedef std::shared_ptr<const PointSource>   ConstPtr;

    PointSource(const Position &position);
    PointSource(const Position &position, const Stokes &stokes);

    virtual const Position &position() const;
    void setPosition(const Position &position);

    void setStokes(const Stokes &stokes);

    template <typename T>
    void setSpectralTerms(double refFreq, bool isLogarithmic, T first, T last);

    void setRotationMeasure(double fraction, double angle, double rm);

    Stokes stokes(double freq) const;

    virtual void accept(ModelComponentVisitor &visitor) const;

private:
    bool hasSpectralTerms() const;
    bool hasRotationMeasure() const;

    Position        itsPosition;
    Stokes          itsStokes;
    double          itsRefFreq;
    std::vector<double>  itsSpectralTerms;
    double          itsPolarizedFraction;
    double          itsPolarizationAngle;
    double          itsRotationMeasure;
    bool            itsHasRotationMeasure;
    bool            itsHasLogarithmicSI;
};

// @}

// -------------------------------------------------------------------------- //
// - Implementation: PointSource                                            - //
// -------------------------------------------------------------------------- //

template <typename T>
void PointSource::setSpectralTerms(double refFreq, bool isLogarithmic, T first, T last)
{
    itsRefFreq = refFreq;
    itsHasLogarithmicSI = isLogarithmic;
    itsSpectralTerms.clear();
    itsSpectralTerms.insert(itsSpectralTerms.begin(), first, last);
}

inline const Position &PointSource::position() const
{
    return itsPosition;
}

} //# namespace DPPP
} //# namespace LOFAR

#endif
