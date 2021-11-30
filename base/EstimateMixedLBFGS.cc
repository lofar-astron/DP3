// EstimateMixed.cc: Estimate Jones matrices for several directions
// simultaneously. A separate data stream is used for each direction. The
// mixing coefficients quantify the influence of each direction on each of the
// other directions (including time and frequency smearing). The solving is done
// using the LBFGS algorithm, with a robust noise model.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "EstimateMixed.h"

#include <casacore/scimath/Fitting/LSQFit.h>

#include "../common/StreamUtil.h"  ///

#include <iostream>
#ifdef HAVE_LIBDIRAC
#include <Dirac.h>
#endif /* HAVE_LIBDIRAC */

using dp3::common::operator<<;

namespace dp3 {
namespace base {

namespace {
#ifdef HAVE_LIBDIRAC
struct LBFGSData {
  std::size_t n_direction;
  std::size_t n_station;
  std::size_t n_baseline;
  std::size_t n_channel;
  const_cursor<Baseline> baselines;
  std::vector<const_cursor<fcomplex>> data;
  std::vector<const_cursor<dcomplex>> model;
  const_cursor<bool> flag;
  const_cursor<float> weight;
  const_cursor<dcomplex> mix;
  double robust_nu;
  LBFGSData(std::size_t n_dir_, std::size_t n_st_, std::size_t n_base_,
            std::size_t n_chan_, const_cursor<Baseline> baselines_,
            std::vector<const_cursor<fcomplex>> data_,
            std::vector<const_cursor<dcomplex>> model_,
            const_cursor<bool> flag_, const_cursor<float> weight_,
            const_cursor<dcomplex> mix_, double robust_nu_ = 2.0)
      : n_direction(n_dir_),
        n_station(n_st_),
        n_baseline(n_base_),
        n_channel(n_chan_),
        baselines(baselines_),
        data(data_),
        model(model_),
        flag(flag_),
        weight(weight_),
        mix(mix_),
        robust_nu(robust_nu_){};
};

// cost function
// unknowns: mx1 vector
double cost_func(double *unknowns, int m, void *adata) {
  assert(adata);
  LBFGSData *t = (LBFGSData *)adata;
  assert(t && t->data.size() == t->n_direction &&
         t->model.size() == t->n_direction);

  const std::size_t n_unknowns = t->n_direction * t->n_station * 4 * 2;
  assert(static_cast<std::size_t>(m) == n_unknowns);

  // Allocate space for intermediate results.
  std::vector<dcomplex> M(t->n_direction * 4);

  double fcost = 0.0;

  // Iterate over baselines
  for (std::size_t bl = 0; bl < t->n_baseline; ++bl) {
    const std::size_t p = t->baselines->first;
    const std::size_t q = t->baselines->second;

    if (p != q) {
      for (std::size_t ch = 0; ch < t->n_channel; ++ch) {
        for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
          // Jones matrix for station P.
          const double *Jp = &(unknowns[dr * t->n_station * 8 + p * 8]);
          const dcomplex Jp_00(Jp[0], Jp[1]);
          const dcomplex Jp_01(Jp[2], Jp[3]);
          const dcomplex Jp_10(Jp[4], Jp[5]);
          const dcomplex Jp_11(Jp[6], Jp[7]);

          // Jones matrix for station Q, conjugated.
          const double *Jq = &(unknowns[dr * t->n_station * 8 + q * 8]);
          const dcomplex Jq_00(Jq[0], -Jq[1]);
          const dcomplex Jq_01(Jq[2], -Jq[3]);
          const dcomplex Jq_10(Jq[4], -Jq[5]);
          const dcomplex Jq_11(Jq[6], -Jq[7]);

          // Fetch model visibilities for the current direction.
          const dcomplex xx = t->model[dr][0];
          const dcomplex xy = t->model[dr][1];
          const dcomplex yx = t->model[dr][2];
          const dcomplex yy = t->model[dr][3];

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
          M[dr * 4] = Jp_00 * Jq_00xx_01xy + Jp_01 * Jq_00yx_01yy;
          M[dr * 4 + 1] = Jp_00 * Jq_10xx_11xy + Jp_01 * Jq_10yx_11yy;
          M[dr * 4 + 2] = Jp_10 * Jq_00xx_01xy + Jp_11 * Jq_00yx_01yy;
          M[dr * 4 + 3] = Jp_10 * Jq_10xx_11xy + Jp_11 * Jq_10yx_11yy;
        }

        for (std::size_t cr = 0; cr < 4; ++cr)  // correlation: 00,01,10,11
        {
          if (!(t->flag[cr])) {
            const double mwt = static_cast<double>(t->weight[cr]);
            for (std::size_t tg = 0; tg < t->n_direction; ++tg) {
              dcomplex visibility(0.0, 0.0);
              for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
                // Look-up mixing weight.
                const dcomplex mix_weight = *(t->mix);

                // Weight model visibility.
                visibility += mix_weight * M[dr * 4 + cr];

                // Move to next source direction.
                t->mix.forward(1);
              }  // Source directions.

              // Compute the residual.
              dcomplex residual =
                  static_cast<dcomplex>(t->data[tg][cr]) - visibility;

              // sum up cost
              // For reference: Gaussian cost is
              // fcost += mwt * (real(residual) * real(residual) +
              //                  imag(residual) * imag(residual));
              // Robust cost function
              fcost += std::log(1.0 + mwt * real(residual) * real(residual) /
                                          t->robust_nu);
              fcost += std::log(1.0 + mwt * imag(residual) * imag(residual) /
                                          t->robust_nu);

              // Move to next target direction.
              t->mix.backward(1, t->n_direction);
              t->mix.forward(0);
            }  // Target directions.

            // Reset cursor to the start of the correlation.
            t->mix.backward(0, t->n_direction);
          }

          // Move to the next correlation.
          t->mix.forward(2);
        }  // Correlations.

        // Move to the next channel.
        t->mix.backward(2, 4);
        t->mix.forward(3);

        for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
          t->model[dr].forward(1);
          t->data[dr].forward(1);
        }
        t->flag.forward(1);
        t->weight.forward(1);
      }  // Channels.

      // Reset cursors to the start of the baseline.
      for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
        t->model[dr].backward(1, t->n_channel);
        t->data[dr].backward(1, t->n_channel);
      }
      t->flag.backward(1, t->n_channel);
      t->weight.backward(1, t->n_channel);
      t->mix.backward(3, t->n_channel);
    }

    // Move cursors to the next baseline.
    for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
      t->model[dr].forward(2);
      t->data[dr].forward(2);
    }
    t->flag.forward(2);
    t->weight.forward(2);
    t->mix.forward(4);
    ++(t->baselines);
  }  // Baselines.

  // Reset all cursors for the next iteration.
  for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
    t->model[dr].backward(2, t->n_baseline);
    t->data[dr].backward(2, t->n_baseline);
  }
  t->flag.backward(2, t->n_baseline);
  t->weight.backward(2, t->n_baseline);
  t->mix.backward(4, t->n_baseline);
  t->baselines -= t->n_baseline;

  return fcost;
}
// gradient function
// unknowns: mx1 parameters, grad: mx1 gradient
/* derivative of baseline p-q:
 *   = 2 real{ vec^H(residual_this_baseline)
            * vec(-J_{pm}C_{pqm} J_{qm}^H) }
        where m: chosen direction
        J_{pm},J_{qm} Jones matrices for baseline p-q
        depending on the parameter, J_p or J_q ==> E
        E: zero matrix, except 1 at location of parameter

 In terms of scalars (used below because of the 'mixing' weight)
 each complex scalar term in the cost: ||z||^2 = z^H z where z=[real j.imag]^T
 of the residual therefore, grad of this term = z^H grad(z) = real grad(real) +
 imag grad(imag) Residual = data- mixing * (J C J^H), expressed for each of the
 4 scalars (2x2 mat)
*/
void grad_func(double *unknowns, double *grad, int m, void *adata) {
  assert(adata);
  LBFGSData *t = (LBFGSData *)adata;
  assert(t && t->data.size() == t->n_direction &&
         t->model.size() == t->n_direction);

  const std::size_t n_unknowns = t->n_direction * t->n_station * 4 * 2;
  assert(static_cast<std::size_t>(m) == n_unknowns);
  const std::size_t n_partial =
      t->n_direction * 8;  // note: this for real,imag separately

  // Allocate space for intermediate results.
  std::vector<dcomplex> M(t->n_direction * 4);
  std::vector<dcomplex> dM(t->n_direction * 16);
  std::vector<double> dR(n_partial);
  std::vector<double> dI(n_partial);
  std::vector<unsigned int> dIndex(4 * n_partial);  // 4 correlations
  // set grad to zero
  std::fill_n(grad, n_unknowns, 0.0);

  // Iterate over baselines and accumulate gradient
  for (std::size_t bl = 0; bl < t->n_baseline; ++bl) {
    const std::size_t p = t->baselines->first;
    const std::size_t q = t->baselines->second;

    if (p != q) {
      // Create partial derivative index for current baseline.
      makeIndex(t->n_direction, t->n_station, *(t->baselines), &(dIndex[0]));

      for (std::size_t ch = 0; ch < t->n_channel; ++ch) {
        for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
          // Jones matrix for station P.
          const double *Jp = &(unknowns[dr * t->n_station * 8 + p * 8]);
          const dcomplex Jp_00(Jp[0], Jp[1]);
          const dcomplex Jp_01(Jp[2], Jp[3]);
          const dcomplex Jp_10(Jp[4], Jp[5]);
          const dcomplex Jp_11(Jp[6], Jp[7]);

          // Jones matrix for station Q, conjugated.
          const double *Jq = &(unknowns[dr * t->n_station * 8 + q * 8]);
          const dcomplex Jq_00(Jq[0], -Jq[1]);
          const dcomplex Jq_01(Jq[2], -Jq[3]);
          const dcomplex Jq_10(Jq[4], -Jq[5]);
          const dcomplex Jq_11(Jq[6], -Jq[7]);

          // Fetch model coherencies for the current direction.
          const dcomplex xx = t->model[dr][0];
          const dcomplex xy = t->model[dr][1];
          const dcomplex yx = t->model[dr][2];
          const dcomplex yy = t->model[dr][3];

          // Precompute terms involving conj(Jq) and the model
          // visibilities.
          const dcomplex Jq_00xx_01xy = Jq_00 * xx + Jq_01 * xy;
          const dcomplex Jq_00yx_01yy = Jq_00 * yx + Jq_01 * yy;
          const dcomplex Jq_10xx_11xy = Jq_10 * xx + Jq_11 * xy;
          const dcomplex Jq_10yx_11yy = Jq_10 * yx + Jq_11 * yy;

          // Precompute (Jp x conj(Jq)) * vec(model), where 'x'
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

        for (std::size_t cr = 0; cr < 4; ++cr)  // correlation: 00,01,10,11
        {
          if (!(t->flag[cr])) {
            const double mwt = static_cast<double>(t->weight[cr]);
            for (std::size_t tg = 0; tg < t->n_direction; ++tg) {
              dcomplex visibility(0.0, 0.0);
              for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
                // Look-up mixing weight.
                const dcomplex mix_weight = *(t->mix);

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
                t->mix.forward(1);
              }  // Source directions.

              // Compute the residual.
              dcomplex residual =
                  static_cast<dcomplex>(t->data[tg][cr]) - visibility;

              // accumulate gradient (for this correlation 'cr')
              // For reference, gradient for Gaussian noise is:
              //  grad[dIndex[cr * n_partial + ci]] +=
              //      2.0 * dR[ci] * mwt * real(residual);
              //  grad[dIndex[cr * n_partial + ci]] +=
              //      2.0 * dI[ci] * mwt * imag(residual);
              for (std::size_t ci = 0; ci < n_partial; ci++) {
                grad[dIndex[cr * n_partial + ci]] -=
                    2.0 * dR[ci] * mwt * real(residual) /
                    (t->robust_nu + mwt * real(residual) * real(residual));
                grad[dIndex[cr * n_partial + ci]] -=
                    2.0 * dI[ci] * mwt * imag(residual) /
                    (t->robust_nu + mwt * imag(residual) * imag(residual));
              }

              // Move to next target direction.
              t->mix.backward(1, t->n_direction);
              t->mix.forward(0);
            }  // Target directions.

            // Reset cursor to the start of the correlation.
            t->mix.backward(0, t->n_direction);
          }

          // Move to the next correlation.
          t->mix.forward(2);
        }  // Correlations.

        // Move to the next channel.
        t->mix.backward(2, 4);
        t->mix.forward(3);

        for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
          t->model[dr].forward(1);
          t->data[dr].forward(1);
        }
        t->flag.forward(1);
        t->weight.forward(1);
      }  // Channels.

      // Reset cursors to the start of the baseline.
      for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
        t->model[dr].backward(1, t->n_channel);
        t->data[dr].backward(1, t->n_channel);
      }
      t->flag.backward(1, t->n_channel);
      t->weight.backward(1, t->n_channel);
      t->mix.backward(3, t->n_channel);
    }

    // Move cursors to the next baseline.
    for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
      t->model[dr].forward(2);
      t->data[dr].forward(2);
    }
    t->flag.forward(2);
    t->weight.forward(2);
    t->mix.forward(4);
    ++(t->baselines);
  }  // Baselines.

  // Reset all cursors for the next iteration.
  for (std::size_t dr = 0; dr < t->n_direction; ++dr) {
    t->model[dr].backward(2, t->n_baseline);
    t->data[dr].backward(2, t->n_baseline);
  }
  t->flag.backward(2, t->n_baseline);
  t->weight.backward(2, t->n_baseline);
  t->mix.backward(4, t->n_baseline);
  t->baselines -= t->n_baseline;
}
#endif /* HAVE_LIBDIRAC */
}  // Unnamed namespace.

#ifdef HAVE_LIBDIRAC
bool estimate(std::size_t n_direction, std::size_t n_station,
              std::size_t n_baseline, std::size_t n_channel,
              const_cursor<Baseline> baselines,
              std::vector<const_cursor<fcomplex>> data,
              std::vector<const_cursor<dcomplex>> model,
              const_cursor<bool> flag, const_cursor<float> weight,
              const_cursor<dcomplex> mix, double *unknowns,
              std::size_t lbfgs_mem, double robust_nu, std::size_t max_iter) {
  LBFGSData t(n_direction, n_station, n_baseline, n_channel, baselines, data,
              model, flag, weight, mix, robust_nu);

  /* for full batch mode, last argument is NULL */
  /* LBFGS memory size lbfgs_mem */
  const int retval =
      lbfgs_fit(cost_func, grad_func, unknowns, n_direction * n_station * 4 * 2,
                max_iter, lbfgs_mem, (void *)&t, nullptr);

  return retval == 0;
}
#endif /* HAVE_LIBDIRAC */

}  // namespace base
}  // namespace dp3
