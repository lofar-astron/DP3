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
    for(size_t bl = 0; bl < nBaseline; ++bl)
    {
        const size_t p = baselines->first;
        const size_t q = baselines->second;

        if(p != q)
        {
            // Jones matrix for station P.
            coeff.forward(1, p);
            const dcomplex Jp_00(coeff[0], coeff[1]);
            const dcomplex Jp_01(coeff[2], coeff[3]);
            const dcomplex Jp_10(coeff[4], coeff[5]);
            const dcomplex Jp_11(coeff[6], coeff[7]);
            coeff.backward(1, p);

            // Jones matrix for station Q, conjugated.
            coeff.forward(1, q);
            const dcomplex Jq_00(coeff[0], -coeff[1]);
            const dcomplex Jq_01(coeff[2], -coeff[3]);
            const dcomplex Jq_10(coeff[4], -coeff[5]);
            const dcomplex Jq_11(coeff[6], -coeff[7]);
            coeff.backward(1, q);

            // Compute (Jp x conj(Jq)) * vec(data), where 'x' denotes the
            // Kronecker product.
            for(size_t ch = 0; ch < nChannel; ++ch)
            {
                // Fetch visibilities.
                const dcomplex xx = data[0];
                const dcomplex xy = data[1];
                const dcomplex yx = data[2];
                const dcomplex yy = data[3];

                // Precompute terms involving conj(Jq) and data. Each term
                // appears twice in the computation of (Jp x conj(Jq))
                // * vec(data).
                const dcomplex Jq_00xx_01xy = Jq_00 * xx + Jq_01 * xy;
                const dcomplex Jq_00yx_01yy = Jq_00 * yx + Jq_01 * yy;
                const dcomplex Jq_10xx_11xy = Jq_10 * xx + Jq_11 * xy;
                const dcomplex Jq_10yx_11yy = Jq_10 * yx + Jq_11 * yy;

                // Compute (Jp x conj(Jq)) * vec(data) from the precomputed
                // terms.
                data[0] = Jp_00 * Jq_00xx_01xy + Jp_01 * Jq_00yx_01yy;
                data[1] = Jp_00 * Jq_10xx_11xy + Jp_01 * Jq_10yx_11yy;
                data[2] = Jp_10 * Jq_00xx_01xy + Jp_11 * Jq_00yx_01yy;
                data[3] = Jp_10 * Jq_10xx_11xy + Jp_11 * Jq_10yx_11yy;

                // Move to the next channel.
                data.forward(1);
            } // Channels.

            // Reset cursor to the beginning of the current baseline.
            data.backward(1, nChannel);
        }

        // Move to the next baseline.
        data.forward(2);
        ++baselines;
    } // Baselines.
}

} //# namespace DPPP
} //# namespace LOFAR
