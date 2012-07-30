//# Patch.cc: A set of sources for which direction dependent effects are assumed to be equal.
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
#include <DPPP/Patch.h>
#include <DPPP/ModelComponentVisitor.h>
#include <Common/lofar_math.h>

namespace LOFAR
{
namespace DPPP
{

Patch::const_iterator Patch::begin() const
{
    return itsComponents.begin();
}

Patch::const_iterator Patch::end() const
{
    return itsComponents.end();
}

void Patch::computePosition()
{
    itsPosition = Position();

    if(!itsComponents.empty())
    {
        double x = 0.0, y = 0.0, z = 0.0;
        for(const_iterator it = begin(), it_end = end(); it != it_end; ++it)
        {
            const Position &position = (*it)->position();
            double cosDec = cos(position[1]);
            x += cos(position[0]) * cosDec;
            y += sin(position[0]) * cosDec;
            z += sin(position[1]);
        }

        x /= itsComponents.size();
        y /= itsComponents.size();
        z /= itsComponents.size();

        itsPosition[0] = atan2(y, x);
        itsPosition[1] = asin(z);
    }
}

} //# namespace DPPP
} //# namespace LOFAR
