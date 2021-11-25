// EstimateMixed.cc: Estimate Jones matrices for several directions
// simultaneously. A separate data stream is used for each direction. The
// mixing coefficients quantify the influence of each direction on each of the
// other directions (including time and frequency smearing).
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "EstimateMixed.h"

#include <casacore/scimath/Fitting/LSQFit.h>

#include "../common/StreamUtil.h"  ///

#include <iostream>
using dp3::common::operator<<;

namespace dp3 {
namespace base {

// Compute a map that contains the index of the unknowns related to the
// specified baseline in the list of all unknowns.
void makeIndex(std::size_t n_direction, std::size_t n_station,
               const Baseline &baseline, unsigned int *index) {
  const std::size_t n_correlation = 4;
  for (std::size_t cr = 0; cr < n_correlation; ++cr) {
    std::size_t idx0 = baseline.first * 8 + (cr / 2) * 4;  // row of P
    std::size_t idx1 =
        baseline.second * 8 + (cr % 2) * 4;  // column of Q (row of Q^H)

    for (std::size_t dr = 0; dr < n_direction; ++dr) {
      *index++ = idx0;
      *index++ = idx0 + 1;
      *index++ = idx0 + 2;
      *index++ = idx0 + 3;

      *index++ = idx1;
      *index++ = idx1 + 1;
      *index++ = idx1 + 2;
      *index++ = idx1 + 3;

      idx0 += n_station * 8;
      idx1 += n_station * 8;
    }
  }
}

bool estimate(size_t nDirection, size_t nStation, size_t nBaseline,
              size_t nChannel, const_cursor<Baseline> baselines,
              vector<const_cursor<fcomplex>> data,
              vector<const_cursor<dcomplex>> model, const_cursor<bool> flag,
              const_cursor<float> weight, const_cursor<dcomplex> mix,
              double *unknowns, size_t maxiter) {
  assert(data.size() == nDirection && model.size() == nDirection);

  // Initialize LSQ solver.
  const size_t nUnknowns = nDirection * nStation * 4 * 2;
  casacore::LSQFit solver(nUnknowns);

  // Each visibility provides information about two (complex) unknowns per
  // station per direction. A visibility is measured by a specific
  // interferometer, which is the combination of two stations. Thus, in total
  // each visibility provides information about (no. of directions) x 2 x 2
  // x 2 (scalar) unknowns = (no. of directions) x 8. For each of these
  // unknowns the value of the partial derivative of the model with respect
  // to the unknown has to be computed.
  const size_t nPartial = nDirection * 8;
  vector<unsigned int> dIndex(4 * nPartial);

  // Allocate space for intermediate results.
  vector<dcomplex> M(nDirection * 4), dM(nDirection * 16);
  vector<double> dR(nPartial), dI(nPartial);

  // Iterate until convergence.
  size_t nIterations = 0;
  while (!solver.isReady() && nIterations < maxiter) {
    for (size_t bl = 0; bl < nBaseline; ++bl) {
      const size_t p = baselines->first;
      const size_t q = baselines->second;

      if (p != q) {
        // Create partial derivative index for current baseline.
        makeIndex(nDirection, nStation, *baselines, &(dIndex[0]));

        for (size_t ch = 0; ch < nChannel; ++ch) {
          for (size_t dr = 0; dr < nDirection; ++dr) {
            // Jones matrix for station P.
            const double *Jp = &(unknowns[dr * nStation * 8 + p * 8]);
            const dcomplex Jp_00(Jp[0], Jp[1]);
            const dcomplex Jp_01(Jp[2], Jp[3]);
            const dcomplex Jp_10(Jp[4], Jp[5]);
            const dcomplex Jp_11(Jp[6], Jp[7]);

            // Jones matrix for station Q, conjugated.
            const double *Jq = &(unknowns[dr * nStation * 8 + q * 8]);
            const dcomplex Jq_00(Jq[0], -Jq[1]);
            const dcomplex Jq_01(Jq[2], -Jq[3]);
            const dcomplex Jq_10(Jq[4], -Jq[5]);
            const dcomplex Jq_11(Jq[6], -Jq[7]);

            // Fetch model visibilities for the current direction.
            const dcomplex xx = model[dr][0];
            const dcomplex xy = model[dr][1];
            const dcomplex yx = model[dr][2];
            const dcomplex yy = model[dr][3];

            // Precompute terms involving conj(Jq) and the model
            // visibilities.
            const dcomplex Jq_00xx_01xy = Jq_00 * xx + Jq_01 * xy;
            const dcomplex Jq_00yx_01yy = Jq_00 * yx + Jq_01 * yy;
            const dcomplex Jq_10xx_11xy = Jq_10 * xx + Jq_11 * xy;
            const dcomplex Jq_10yx_11yy = Jq_10 * yx + Jq_11 * yy;

            // Precompute (Jp x conj(Jq)) * vec(data), where 'x'
            // denotes the Kronecker product. This is the model
            // visibility for the current direction, with the
            // current Jones matrix estimates applied. This is
            // stored in M.
            // Also, precompute the partial derivatives of M with
            // respect to all 16 parameters (i.e. 2 Jones matrices
            // Jp and Jq, 4 complex scalars per Jones matrix, 2 real
            // scalars per complex scalar, 2 * 4 * 2 = 16). These
            // partial derivatives are stored in dM.
            M[dr * 4] = Jp_00 * Jq_00xx_01xy + Jp_01 * Jq_00yx_01yy;
            dM[dr * 16] = Jq_00xx_01xy;                 // dM_00/dJp_00
            dM[dr * 16 + 1] = Jq_00yx_01yy;             // dM_00/dJp_01
            dM[dr * 16 + 2] = Jp_00 * xx + Jp_01 * yx;  // dM_00/dJq_00
            dM[dr * 16 + 3] = Jp_00 * xy + Jp_01 * yy;  // dM_00/dJq_01

            M[dr * 4 + 1] = Jp_00 * Jq_10xx_11xy + Jp_01 * Jq_10yx_11yy;
            dM[dr * 16 + 4] = Jq_10xx_11xy;     // dM_01/dJp_00
            dM[dr * 16 + 5] = Jq_10yx_11yy;     // dM_01/dJp_01
            dM[dr * 16 + 6] = dM[dr * 16 + 2];  // dM_01/dJq_10
            dM[dr * 16 + 7] = dM[dr * 16 + 3];  // dM_01/dJq_11

            M[dr * 4 + 2] = Jp_10 * Jq_00xx_01xy + Jp_11 * Jq_00yx_01yy;
            dM[dr * 16 + 8] = dM[dr * 16];               // dM_10/dJp_10
            dM[dr * 16 + 9] = dM[dr * 16 + 1];           // dM_10/dJp_11
            dM[dr * 16 + 10] = Jp_10 * xx + Jp_11 * yx;  // dM_10/dJq_00
            dM[dr * 16 + 11] = Jp_10 * xy + Jp_11 * yy;  // dM_10/dJq_01

            M[dr * 4 + 3] = Jp_10 * Jq_10xx_11xy + Jp_11 * Jq_10yx_11yy;
            dM[dr * 16 + 12] = dM[dr * 16 + 4];   // dM_11/dJp_10
            dM[dr * 16 + 13] = dM[dr * 16 + 5];   // dM_11/dJp_11
            dM[dr * 16 + 14] = dM[dr * 16 + 10];  // dM_11/dJq_10
            dM[dr * 16 + 15] = dM[dr * 16 + 11];  // dM_11/dJq_11
          }

          for (size_t cr = 0; cr < 4; ++cr)  // correlation: 00,01,10,11
          {
            if (!flag[cr]) {
              for (size_t tg = 0; tg < nDirection; ++tg) {
                dcomplex visibility(0.0, 0.0);
                for (size_t dr = 0; dr < nDirection; ++dr) {
                  // Look-up mixing weight.
                  const dcomplex mix_weight = *mix;

                  // Weight model visibility.
                  visibility += mix_weight * M[dr * 4 + cr];

                  // Compute weighted partial derivatives.
                  dcomplex derivative(0.0, 0.0);
                  derivative = mix_weight * dM[dr * 16 + cr * 4];
                  dR[dr * 8] = real(derivative);  // for cr==0: Re(d/dRe(p_00)))
                  dI[dr * 8] = imag(derivative);  // for cr==0: Re(d/dIm(p_00)))
                  dR[dr * 8 + 1] =
                      -imag(derivative);  // for cr==0: Im(d/dRe(p_00)))
                  dI[dr * 8 + 1] =
                      real(derivative);  // for cr==0: Im(d/dIm(p_00)))

                  derivative = mix_weight * dM[dr * 16 + cr * 4 + 1];
                  dR[dr * 8 + 2] =
                      real(derivative);  // for cr==0: Re(d/dRe(p_01)))
                  dI[dr * 8 + 2] =
                      imag(derivative);  // for cr==0: Re(d/dIm(p_01)))
                  dR[dr * 8 + 3] =
                      -imag(derivative);  // for cr==0: Im(d/dRe(p_01)))
                  dI[dr * 8 + 3] =
                      real(derivative);  // for cr==0: Im(d/dIm(p_01)))

                  derivative = mix_weight * dM[dr * 16 + cr * 4 + 2];
                  dR[dr * 8 + 4] =
                      real(derivative);  // for cr==0: Re(d/dRe(q_00)))
                  dI[dr * 8 + 4] =
                      imag(derivative);  // for cr==0: Re(d/dIm(q_00)))
                  dR[dr * 8 + 5] =
                      imag(derivative);  // for cr==0: Im(d/dRe(q_00)))
                  dI[dr * 8 + 5] =
                      -real(derivative);  // for cr==0: Im(d/dIm(q_00)))

                  derivative = mix_weight * dM[dr * 16 + cr * 4 + 3];
                  dR[dr * 8 + 6] =
                      real(derivative);  // for cr==0: Re(d/dRe(q_01)))
                  dI[dr * 8 + 6] =
                      imag(derivative);  // for cr==0: Re(d/dIm(q_01)))
                  dR[dr * 8 + 7] =
                      imag(derivative);  // for cr==0: Im(d/dRe(q_01)))
                  dI[dr * 8 + 7] =
                      -real(derivative);  // for cr==0: Im(d/dIm(q_01)))

                  // Move to next source direction.
                  mix.forward(1);
                }  // Source directions.

                // Compute the residual.
                dcomplex residual =
                    static_cast<dcomplex>(data[tg][cr]) - visibility;

                // Update the normal equations.
                solver.makeNorm(nPartial, &(dIndex[cr * nPartial]), &(dR[0]),
                                static_cast<double>(weight[cr]),
                                real(residual));
                solver.makeNorm(nPartial, &(dIndex[cr * nPartial]), &(dI[0]),
                                static_cast<double>(weight[cr]),
                                imag(residual));

                // Move to next target direction.
                mix.backward(1, nDirection);
                mix.forward(0);
              }  // Target directions.

              // Reset cursor to the start of the correlation.
              mix.backward(0, nDirection);
            }

            // Move to the next correlation.
            mix.forward(2);
          }  // Correlations.

          // Move to the next channel.
          mix.backward(2, 4);
          mix.forward(3);

          for (size_t dr = 0; dr < nDirection; ++dr) {
            model[dr].forward(1);
            data[dr].forward(1);
          }
          flag.forward(1);
          weight.forward(1);
        }  // Channels.

        // Reset cursors to the start of the baseline.
        for (size_t dr = 0; dr < nDirection; ++dr) {
          model[dr].backward(1, nChannel);
          data[dr].backward(1, nChannel);
        }
        flag.backward(1, nChannel);
        weight.backward(1, nChannel);
        mix.backward(3, nChannel);
      }

      // Move cursors to the next baseline.
      for (size_t dr = 0; dr < nDirection; ++dr) {
        model[dr].forward(2);
        data[dr].forward(2);
      }
      flag.forward(2);
      weight.forward(2);
      mix.forward(4);
      ++baselines;
    }  // Baselines.

    // Reset all cursors for the next iteration.
    for (size_t dr = 0; dr < nDirection; ++dr) {
      model[dr].backward(2, nBaseline);
      data[dr].backward(2, nBaseline);
    }
    flag.backward(2, nBaseline);
    weight.backward(2, nBaseline);
    mix.backward(4, nBaseline);
    baselines -= nBaseline;

    // Perform LSQ iteration.
    casacore::uInt rank;
    [[maybe_unused]] bool status = solver.solveLoop(rank, unknowns, true);
    assert(status);

    // Update iteration count.
    ++nIterations;
  }

  const bool converged = (solver.isReady() == casacore::LSQFit::SOLINCREMENT ||
                          solver.isReady() == casacore::LSQFit::DERIVLEVEL);
  return converged;
}

}  // namespace base
}  // namespace dp3
