//# SubtractMixed.cc: Subtract visibilities from a buffer after weighting by
//# mixing coefficients.
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
#include <DPPP/SubtractMixed.h>

namespace LOFAR
{
namespace DPPP
{

void subtract(size_t nBaseline, size_t nChannel,
    const_cursor<Baseline> baselines, cursor<fcomplex> data,
    const_cursor<dcomplex> model, const_cursor<dcomplex> weight)
{
    for(size_t bl = 0; bl < nBaseline; ++bl)
    {
        const size_t p = baselines->first;
        const size_t q = baselines->second;

        if(p != q)
        {
            for(size_t ch = 0; ch < nChannel; ++ch)
            {
                // Subtract weighted model from data.
                *data -= static_cast<fcomplex>((*weight) * (*model));
                ++weight;
                ++model;
                ++data;
                *data -= static_cast<fcomplex>((*weight) * (*model));
                ++weight;
                ++model;
                ++data;
                *data -= static_cast<fcomplex>((*weight) * (*model));
                ++weight;
                ++model;
                ++data;
                *data -= static_cast<fcomplex>((*weight) * (*model));
                ++weight;
                ++model;
                ++data;

                // Move to the next channel.
                weight -= 4;
                weight.forward(1);
                model -= 4;
                model.forward(1);
                data -= 4;
                data.forward(1);
            } // Channels.

            weight.backward(1, nChannel);
            model.backward(1, nChannel);
            data.backward(1, nChannel);
        }

        // Move to the next baseline.
        weight.forward(2);
        model.forward(2);
        data.forward(2);
        ++baselines;
    } // Baselines.
}

} //# namespace DPPP
} //# namespace LOFAR
