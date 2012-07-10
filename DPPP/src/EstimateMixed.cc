//# EstimateMixed.cc: Estimate Jones matrices for several directions
//# simultaneously. A separate data stream is used for each direction. The
//# mixing coefficients quantify the influence of each direction on each of the
//# other directions (including time and frequency smearing).
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
#include <DPPP/EstimateMixed.h>
#include <Common/LofarLogger.h>
#include <scimath/Fitting/LSQFit.h>
#include <boost/multi_array.hpp>

namespace LOFAR
{
namespace DPPP
{

bool estimate(size_t nDirection, size_t nStation, size_t nBaseline,
    size_t nChannel, const_cursor<Baseline> baselines,
    vector<const_cursor<fcomplex> > data, vector<const_cursor<dcomplex> > model,
    const_cursor<bool> flag, const_cursor<float> weight,
    const_cursor<dcomplex> mix, double *unknowns, double *errors)
{
    ASSERT(data.size() == nDirection && model.size() == nDirection);

    // Initialize LSQ solver.
    const size_t nUnknowns = nDirection * nStation * 4 * 2;
    casa::LSQFit solver(nUnknowns);

    // Each visibility provides information about two (complex) unknowns per
    // station per direction. A visibility is measured by a specific
    // interferometer, which is the combination of two stations. Thus, in total
    // each visibility provides information about (no. of directions) x 2 x 2
    // x 2 (scalar) unknowns = (no. of directions) x 8. For each of these
    // unknowns the value of the partial derivative of the model with respect
    // to the unknow has to be computed.
    const unsigned int nPartial = nDirection * 8;

    // Construct partial derivative index template for each correlation.
    boost::multi_array<unsigned int, 2> dIndexTemplate(boost::extents[4]
        [nPartial]);
    for(size_t cr = 0; cr < 4; ++cr)
    {
        size_t idx0 = (cr / 2) * 4;
        size_t idx1 = (cr & 1) * 4;

        for(size_t dr = 0; dr < nDirection; ++dr)
        {
            dIndexTemplate[cr][dr * 8 + 0] = idx0 + 0;
            dIndexTemplate[cr][dr * 8 + 1] = idx0 + 1;
            dIndexTemplate[cr][dr * 8 + 2] = idx0 + 2;
            dIndexTemplate[cr][dr * 8 + 3] = idx0 + 3;
            dIndexTemplate[cr][dr * 8 + 4] = idx1 + 0;
            dIndexTemplate[cr][dr * 8 + 5] = idx1 + 1;
            dIndexTemplate[cr][dr * 8 + 6] = idx1 + 2;
            dIndexTemplate[cr][dr * 8 + 7] = idx1 + 3;
            idx0 += nStation * 8;
            idx1 += nStation * 8;
        }
    }

    // Allocate space for intermediate results.
    boost::multi_array<dcomplex, 2> M(boost::extents[nDirection][4]);
    boost::multi_array<dcomplex, 3> dM(boost::extents[nDirection][4][4]);
    boost::multi_array<double, 1> dR(boost::extents[nPartial]);
    boost::multi_array<double, 1> dI(boost::extents[nPartial]);
    boost::multi_array<unsigned int, 2> dIndex(boost::extents[4][nPartial]);

    dcomplex coherence;
    dcomplex Jp_00, Jp_01, Jp_10, Jp_11;
    dcomplex Jq_00, Jq_01, Jq_10, Jq_11;
    dcomplex Jp_00_s0, Jp_10_s0, Jp_00_s1, Jp_10_s1, Jp_01_s2, Jp_11_s2,
        Jp_01_s3, Jp_11_s3;
    dcomplex Jq_00_s0, Jq_10_s0, Jq_01_s1, Jq_11_s1, Jq_00_s2, Jq_10_s2,
        Jq_01_s3, Jq_11_s3;

    // Iterate until convergence.
    size_t nIterations = 0;
    while(!solver.isReady() && nIterations < 50)
    {
        for(size_t bl = 0; bl < nBaseline; ++bl)
        {
            const size_t p = baselines->first;
            const size_t q = baselines->second;

            if(p != q)
            {
                // Update partial derivative index for current baseline.
                const size_t offsetP = p * 8;
                const size_t offsetQ = q * 8;
                for(size_t cr = 0; cr < 4; ++cr)
                {
                    for(size_t dr = 0; dr < nDirection; ++dr)
                    {
                        dIndex[cr][dr * 8 + 0] = dIndexTemplate[cr][dr * 8 + 0]
                            + offsetP;
                        dIndex[cr][dr * 8 + 1] = dIndexTemplate[cr][dr * 8 + 1]
                            + offsetP;
                        dIndex[cr][dr * 8 + 2] = dIndexTemplate[cr][dr * 8 + 2]
                            + offsetP;
                        dIndex[cr][dr * 8 + 3] = dIndexTemplate[cr][dr * 8 + 3]
                            + offsetP;

                        dIndex[cr][dr * 8 + 4] = dIndexTemplate[cr][dr * 8 + 4]
                            + offsetQ;
                        dIndex[cr][dr * 8 + 5] = dIndexTemplate[cr][dr * 8 + 5]
                            + offsetQ;
                        dIndex[cr][dr * 8 + 6] = dIndexTemplate[cr][dr * 8 + 6]
                            + offsetQ;
                        dIndex[cr][dr * 8 + 7] = dIndexTemplate[cr][dr * 8 + 7]
                            + offsetQ;
                    }
                }

                for(size_t ch = 0; ch < nChannel; ++ch)
                {
                    dcomplex coherence;
                    for(size_t dr = 0; dr < nDirection; ++dr)
                    {
                        const double *Jp =
                            &(unknowns[dr * nStation * 8 + p * 8]);
                        Jp_00 = dcomplex(Jp[0], Jp[1]);
                        Jp_01 = dcomplex(Jp[2], Jp[3]);
                        Jp_10 = dcomplex(Jp[4], Jp[5]);
                        Jp_11 = dcomplex(Jp[6], Jp[7]);

                        const double *Jq =
                            &(unknowns[dr * nStation * 8 + q * 8]);
                        Jq_00 = dcomplex(Jq[0], -Jq[1]);
                        Jq_01 = dcomplex(Jq[2], -Jq[3]);
                        Jq_10 = dcomplex(Jq[4], -Jq[5]);
                        Jq_11 = dcomplex(Jq[6], -Jq[7]);

                        coherence = (model[dr][0]);
                        Jp_00_s0 = Jp_00 * coherence;
                        Jp_10_s0 = Jp_10 * coherence;
                        Jq_00_s0 = Jq_00 * coherence;
                        Jq_10_s0 = Jq_10 * coherence;

                        coherence = (model[dr][1]);
                        Jp_00_s1 = Jp_00 * coherence;
                        Jp_10_s1 = Jp_10 * coherence;
                        Jq_01_s1 = Jq_01 * coherence;
                        Jq_11_s1 = Jq_11 * coherence;

                        coherence = (model[dr][2]);
                        Jp_01_s2 = Jp_01 * coherence;
                        Jp_11_s2 = Jp_11 * coherence;
                        Jq_00_s2 = Jq_00 * coherence;
                        Jq_10_s2 = Jq_10 * coherence;

                        coherence = (model[dr][3]);
                        Jp_01_s3 = Jp_01 * coherence;
                        Jp_11_s3 = Jp_11 * coherence;
                        Jq_01_s3 = Jq_01 * coherence;
                        Jq_11_s3 = Jq_11 * coherence;

                        M[dr][0] = Jp_00 * (Jq_00_s0 + Jq_01_s1)
                            + Jp_01 * (Jq_00_s2 + Jq_01_s3);
                        dM[dr][0][0] = Jq_00_s0 + Jq_01_s1;
                        dM[dr][0][1] = Jq_00_s2 + Jq_01_s3;
                        dM[dr][0][2] = Jp_00_s0 + Jp_01_s2;
                        dM[dr][0][3] = Jp_00_s1 + Jp_01_s3;

                        M[dr][1] = Jp_00 * (Jq_10_s0 + Jq_11_s1)
                            + Jp_01 * (Jq_10_s2 + Jq_11_s3);
                        dM[dr][1][0] = Jq_10_s0 + Jq_11_s1;
                        dM[dr][1][1] = Jq_10_s2 + Jq_11_s3;
                        dM[dr][1][2] = dM[dr][0][2];
                        dM[dr][1][3] = dM[dr][0][3];

                        M[dr][2] = Jp_10 * (Jq_00_s0 + Jq_01_s1)
                            + Jp_11 * (Jq_00_s2 + Jq_01_s3);
                        dM[dr][2][0] = dM[dr][0][0];
                        dM[dr][2][1] = dM[dr][0][1];
                        dM[dr][2][2] = Jp_10_s0 + Jp_11_s2;
                        dM[dr][2][3] = Jp_10_s1 + Jp_11_s3;

                        M[dr][3] = Jp_10 * (Jq_10_s0 + Jq_11_s1)
                            + Jp_11 * (Jq_10_s2 + Jq_11_s3);
                        dM[dr][3][0] = dM[dr][1][0];
                        dM[dr][3][1] = dM[dr][1][1];
                        dM[dr][3][2] = dM[dr][2][2];
                        dM[dr][3][3] = dM[dr][2][3];
                    }

                    for(size_t cr = 0; cr < 4; ++cr)
                    {
                        if(!flag[cr])
                        {
                            for(size_t tg = 0; tg < nDirection; ++tg)
                            {
                                dcomplex model = 0.0, partial;
                                for(size_t dr = 0; dr < nDirection; ++dr)
                                {
                                    // Look-up mixing term.
                                    const dcomplex term = *mix;

                                    // Update model visibility.
                                    model += term * M[dr][cr];

                                    // Compute partial derivatives.
                                    partial = term * dM[dr][cr][0];
                                    dR[dr * 8] = real(partial);
                                    dI[dr * 8] = imag(partial);
                                    dR[dr * 8 + 1] = -imag(partial);
                                    dI[dr * 8 + 1] = real(partial);

                                    partial = term * dM[dr][cr][1];
                                    dR[dr * 8 + 2] = real(partial);
                                    dI[dr * 8 + 2] = imag(partial);
                                    dR[dr * 8 + 3] = -imag(partial);
                                    dI[dr * 8 + 3] = real(partial);

                                    partial = term * dM[dr][cr][2];
                                    dR[dr * 8 + 4] = real(partial);
                                    dI[dr * 8 + 4] = imag(partial);
                                    dR[dr * 8 + 5] = imag(partial);
                                    dI[dr * 8 + 5] = -real(partial);

                                    partial = term * dM[dr][cr][3];
                                    dR[dr * 8 + 6] = real(partial);
                                    dI[dr * 8 + 6] = imag(partial);
                                    dR[dr * 8 + 7] = imag(partial);
                                    dI[dr * 8 + 7] = -real(partial);

                                    // Move to next source direction.
                                    mix.forward(1);
                                } // Source directions.

                                // Compute the residual.
                                dcomplex residual =
                                    static_cast<dcomplex>(data[tg][cr]) - model;

                                // Update the normal equations.
                                solver.makeNorm(nPartial, &(dIndex[cr][0]),
                                    &(dR[0]), static_cast<double>(weight[cr]),
                                    real(residual));
                                solver.makeNorm(nPartial, &(dIndex[cr][0]),
                                    &(dI[0]), static_cast<double>(weight[cr]),
                                    imag(residual));

                                // Move to next target direction.
                                mix.backward(1, nDirection);
                                mix.forward(0);
                            } // Target directions.

                            // Reset cursor to the start of the correlation.
                            mix.backward(0, nDirection);
                        }

                        // Move to the next correlation.
                        mix.forward(2);
                    } // Correlations.

                    // Move to the next channel.
                    mix.backward(2, 4);
                    mix.forward(3);

                    for(size_t dr = 0; dr < nDirection; ++dr)
                    {
                        model[dr].forward(1);
                        data[dr].forward(1);
                    }
                    flag.forward(1);
                    weight.forward(1);
                } // Channels.

                // Reset cursors to the start of the baseline.
                for(size_t dr = 0; dr < nDirection; ++dr)
                {
                    model[dr].backward(1, nChannel);
                    data[dr].backward(1, nChannel);
                }
                flag.backward(1, nChannel);
                weight.backward(1, nChannel);
                mix.backward(3, nChannel);
            }

            // Move cursors to the next baseline.
            for(size_t dr = 0; dr < nDirection; ++dr)
            {
                model[dr].forward(2);
                data[dr].forward(2);
            }
            flag.forward(2);
            weight.forward(2);
            mix.forward(4);
            ++baselines;
        } // Baselines.

        // Reset all cursors for the next iteration.
        for(size_t dr = 0; dr < nDirection; ++dr)
        {
            model[dr].backward(2, nBaseline);
            data[dr].backward(2, nBaseline);
        }
        flag.backward(2, nBaseline);
        weight.backward(2, nBaseline);
        mix.backward(4, nBaseline);
        baselines -= nBaseline;

        // Perform LSQ iteration.
        casa::uInt rank;
        bool status = solver.solveLoop(rank, unknowns, true);
        ASSERT(status);

        // Update iteration count.
        ++nIterations;
    }

    // Get the estimated error for each unknown from the solver.
    if(errors)
    {
//        boost::multi_array<double, 1> errors_packed(boost::extents[nUnknowns * nUnknowns]);
//        // TODO: Somehow this only returns zero??
//        bool status = solver.getErrors(&(errors_packed[0]));

        vector<double> cov(nUnknowns * nUnknowns);
        bool status = solver.getCovariance(&(cov[0]));
        ASSERT(status);
        for(size_t i = 0; i < nUnknowns; ++i)
        {
            errors[i] = sqrt(abs(cov[i * nUnknowns + i])) * solver.getSD();
        }
    }

    const bool converged = (solver.isReady() == casa::LSQFit::SOLINCREMENT
        || solver.isReady() == casa::LSQFit::DERIVLEVEL);
    return converged;
}

} //# namespace DPPP
} //# namespace LOFAR
