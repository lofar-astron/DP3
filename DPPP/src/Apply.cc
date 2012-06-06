//# Apply.cc: Apply station Jones matrices to a set of visibilities.
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
#include <DPPP/Apply.h>

namespace LOFAR
{
namespace DPPP
{

void apply(size_t nBaseline, size_t nChannel, const_cursor<Baseline> baselines,
    const_cursor<double> coeff, cursor<dcomplex> data)
{
    dcomplex Jp_00, Jp_01, Jp_10, Jp_11;
    dcomplex Jq_00, Jq_01, Jq_10, Jq_11;
    dcomplex Jq_00_s0, Jq_10_s0, Jq_01_s1, Jq_11_s1, Jq_00_s2, Jq_10_s2,
        Jq_01_s3, Jq_11_s3;

    for(size_t bl = 0; bl < nBaseline; ++bl)
    {
        const size_t p = baselines->first;
        const size_t q = baselines->second;

        if(p != q)
        {
            coeff.forward(1, p);
            Jp_00 = dcomplex(coeff[0], coeff[1]);
            Jp_01 = dcomplex(coeff[2], coeff[3]);
            Jp_10 = dcomplex(coeff[4], coeff[5]);
            Jp_11 = dcomplex(coeff[6], coeff[7]);
            coeff.backward(1, p);

            coeff.forward(1, q);
            Jq_00 = dcomplex(coeff[0], -coeff[1]);
            Jq_01 = dcomplex(coeff[2], -coeff[3]);
            Jq_10 = dcomplex(coeff[4], -coeff[5]);
            Jq_11 = dcomplex(coeff[6], -coeff[7]);
            coeff.backward(1, q);

            for(size_t ch = 0; ch < nChannel; ++ch)
            {
                Jq_00_s0 = Jq_00 * (*data);
                Jq_10_s0 = Jq_10 * (*data);
                ++data;

                Jq_01_s1 = Jq_01 * (*data);
                Jq_11_s1 = Jq_11 * (*data);
                ++data;

                Jq_00_s2 = Jq_00 * (*data);
                Jq_10_s2 = Jq_10 * (*data);
                ++data;

                Jq_01_s3 = Jq_01 * (*data);
                Jq_11_s3 = Jq_11 * (*data);
                ++data;
                data -= 4;

                *data = Jp_00 * (Jq_00_s0 + Jq_01_s1)
                    + Jp_01 * (Jq_00_s2 + Jq_01_s3);
                ++data;

                *data = Jp_00 * (Jq_10_s0 + Jq_11_s1)
                    + Jp_01 * (Jq_10_s2 + Jq_11_s3);
                ++data;

                *data = Jp_10 * (Jq_00_s0 + Jq_01_s1)
                    + Jp_11 * (Jq_00_s2 + Jq_01_s3);
                ++data;

                *data = Jp_10 * (Jq_10_s0 + Jq_11_s1)
                    + Jp_11 * (Jq_10_s2 + Jq_11_s3);
                ++data;
                data -= 4;

                // Move to next channel.
                data.forward(1);
            } // Channels.
            data.backward(1, nChannel);
        }

        data.forward(2);
        ++baselines;
    } // Baselines.
}

} //# namespace DPPP
} //# namespace LOFAR
