// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "LBFGSSolver.h"

// Include before Dirac.h to avoid compilation issue with complex's definition.
#include "common/MatrixComplexDouble2x2.h"
#include "common/DiagonalMatrixComplexDouble2x2.h"

#ifdef HAVE_LIBDIRAC
#include <Dirac.h>
// Dirac.h incorrectly defines complex, we undefine it below
// to avoid conflicts with xtensor/xcomplex.hpp
#undef complex
#endif /* HAVE_LIBDIRAC */

#include <aocommon/parallelfor.h>
#include <xtensor/xadapt.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xlayout.hpp>
#include <xtensor/xvectorize.hpp>
#include <xtensor/xview.hpp>

#include <algorithm>
#include <iostream>

namespace dp3 {
namespace ddecal {

#ifdef HAVE_LIBDIRAC
namespace {

struct lbfgs_fulljones_data {
  const SolveData::ChannelBlockData& cb_data;
  std::size_t n_antennas;
  std::size_t n_directions;
  double robust_nu;
  std::size_t start_baseline;
  std::size_t end_baseline;
  lbfgs_fulljones_data(const SolveData::ChannelBlockData& cb_data_,
                       std::size_t n_antennas_, std::size_t n_directions_,
                       double robust_nu_)
      : cb_data(cb_data_),
        n_antennas(n_antennas_),
        n_directions(n_directions_),
        robust_nu(robust_nu_),
        start_baseline(0),
        end_baseline(0){};
};

double FullCost(double* unknowns, int n_uknowns, void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions * 4;
  // Map unknowns into DComplex vector
  // Combine the real and imaginary values that are in the first and second half
  // of the unknowns array, respectively
  std::vector<SolverBase::DComplex> solutions(n_sol);
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  double cost = 0.0;
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const aocommon::MatrixComplexDouble2x2 J1(
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const aocommon::MatrixComplexDouble2x2 J2(
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);
      const aocommon::MatrixComplexDouble2x2 J1_M = J1 * model_data;
      const aocommon::MatrixComplexDouble2x2 J1_M_J2H =
          J1_M * J2.HermTranspose();
      res -= J1_M_J2H;
    }
    // For LS, cost += Norm(res);
    cost += std::log(1.0 + Norm(res) / lbfgs_dat->robust_nu);
  }

  // normalize cost by number of baselines
  return cost / double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline);
}

void FullGradient(double* unknowns, double* gradient, int n_unknowns,
                  void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions * 4;
  std::vector<SolverBase::DComplex> solutions(n_sol);
  // Map unknowns into DComplex vector
  // Combine the real and imaginary values that are in the first and second half
  // of the unknowns array, respectively
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  std::vector<SolverBase::DComplex> complex_gradient(
      n_sol, SolverBase::DComplex(0.0, 0.0));

  const double inv_baselines =
      1.0 / (double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline));
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
    // loop for residual calculation
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const aocommon::MatrixComplexDouble2x2 J1(
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const aocommon::MatrixComplexDouble2x2 J2(
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);
      const aocommon::MatrixComplexDouble2x2 J1_M = J1 * model_data;
      const aocommon::MatrixComplexDouble2x2 J1_M_J2H =
          J1_M * J2.HermTranspose();
      res -= J1_M_J2H;
    }
    // Scale factor for robust gradient, divided by number_of_baselines
    // For LS cost, scale_factor=2.0/number_of_baselines
    const double scale_factor =
        -2.0 * inv_baselines /
        (lbfgs_dat->robust_nu + Norm(res));  //-ve for -ve grad direction
    // loop for grad calculation
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const aocommon::MatrixComplexDouble2x2 J1(
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const aocommon::MatrixComplexDouble2x2 J2(
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);

      const aocommon::MatrixComplexDouble2x2 res_J2 = res * J2;
      const aocommon::MatrixComplexDouble2x2 res_J2_modH =
          res_J2 * model_data.HermTranspose();
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4] +=
          scale_factor * res_J2_modH[0];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 1] +=
          scale_factor * res_J2_modH[1];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 2] +=
          scale_factor * res_J2_modH[2];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 3] +=
          scale_factor * res_J2_modH[3];

      const aocommon::MatrixComplexDouble2x2 resH_J1 = res.HermTranspose() * J1;
      const aocommon::MatrixComplexDouble2x2 resH_J1_mod = resH_J1 * model_data;
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 4] +=
          scale_factor * resH_J1_mod[0];
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 4 + 1] +=
          scale_factor * resH_J1_mod[1];
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 4 + 2] +=
          scale_factor * resH_J1_mod[2];
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 4 + 3] +=
          scale_factor * resH_J1_mod[3];
    }
  }

  // Copy DComplex vector into gradient (real and imaginary parts separately)
  std::transform(complex_gradient.begin(), complex_gradient.end(), gradient,
                 [](const SolverBase::DComplex& z) { return z.real(); });
  std::transform(complex_gradient.begin(), complex_gradient.end(),
                 gradient + n_sol,
                 [](const SolverBase::DComplex& z) { return z.imag(); });
}

double DiagonalCost(double* unknowns, int n_uknowns, void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions * 2;
  std::vector<SolverBase::DComplex> solutions(n_sol);
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  double cost = 0.0;
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const aocommon::DiagonalMatrixComplexDouble2x2 J1(
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 2]);
      const aocommon::DiagonalMatrixComplexDouble2x2 J2(
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 2]);
      const aocommon::MatrixComplexDouble2x2 J1_M = J1 * model_data;
      const aocommon::MatrixComplexDouble2x2 J1_M_J2H =
          J1_M * J2.HermTranspose();
      res -= J1_M_J2H;
    }
    // For LS, cost += Norm(res);
    cost += std::log(1.0 + Norm(res) / lbfgs_dat->robust_nu);
  }

  // normalize cost by number of baselines
  return cost / double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline);
}

void DiagonalGradient(double* unknowns, double* gradient, int n_unknowns,
                      void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions * 2;
  std::vector<SolverBase::DComplex> solutions(n_sol);
  // Map unknowns into DComplex vector
  // Combine the real and imaginary values that are in the first and second half
  // of the unknowns array, respectively
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  std::vector<SolverBase::DComplex> complex_gradient(
      n_sol, SolverBase::DComplex(0.0, 0.0));

  const double inv_baselines =
      1.0 / (double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline));
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
    // loop for residual calculation
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const aocommon::DiagonalMatrixComplexDouble2x2 J1(
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 2]);
      const aocommon::DiagonalMatrixComplexDouble2x2 J2(
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 2]);
      const aocommon::MatrixComplexDouble2x2 J1_M = J1 * model_data;
      const aocommon::MatrixComplexDouble2x2 J1_M_J2H =
          J1_M * J2.HermTranspose();
      res -= J1_M_J2H;
    }
    // Scale factor for robust gradient, divided by number_of_baselines
    // For LS cost, scale_factor=2.0/number_of_baselines
    const double scale_factor =
        -2.0 * inv_baselines /
        (lbfgs_dat->robust_nu + Norm(res));  //-ve for -ve grad direction
    // loop for grad calculation
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      // C
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const SolverBase::DComplex* g_1 =
          &solutions[(antenna1 * lbfgs_dat->n_directions + d) * 2];
      const SolverBase::DComplex* g_2 =
          &solutions[(antenna2 * lbfgs_dat->n_directions + d) * 2];

      // Hadamard product R o C^*
      const aocommon::MatrixComplexDouble2x2 R_Cconj(
          res[0] * std::conj(model_data[0]), res[1] * std::conj(model_data[1]),
          res[2] * std::conj(model_data[2]), res[3] * std::conj(model_data[3]));

      // grad for antenna 2 = - g_1^H (R o C^*)
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 2] +=
          scale_factor *
          (std::conj(g_1[0]) * R_Cconj[0] + std::conj(g_1[1]) * R_Cconj[2]);
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d) * 2 + 1] +=
          scale_factor *
          (std::conj(g_1[0]) * R_Cconj[1] + std::conj(g_1[1]) * R_Cconj[3]);
      // grad for antenna 1 = - g_2^H (R o C^*)^H
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 2] +=
          scale_factor * std::conj(g_2[0] * R_Cconj[0] + g_2[1] * R_Cconj[1]);
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 2 + 1] +=
          scale_factor * std::conj(g_2[0] * R_Cconj[2] + g_2[1] * R_Cconj[3]);
    }
  }

  // Copy DComplex vector into gradient (real and imaginary parts separately)
  std::transform(complex_gradient.begin(), complex_gradient.end(), gradient,
                 [](const SolverBase::DComplex& z) { return z.real(); });
  // Note : return conjugate of gradient
  std::transform(complex_gradient.begin(), complex_gradient.end(),
                 gradient + n_sol,
                 [](const SolverBase::DComplex& z) { return -z.imag(); });
}

double ScalarCost(double* unknowns, int n_uknowns, void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions;
  std::vector<SolverBase::DComplex> solutions(n_sol);
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  double cost = 0.0;
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const SolverBase::DComplex* g_1 =
          &solutions[(antenna1 * lbfgs_dat->n_directions + d)];
      const SolverBase::DComplex* g_2 =
          &solutions[(antenna2 * lbfgs_dat->n_directions + d)];
      res -= model_data * (g_1[0] * std::conj(g_2[0]));
    }
    // For LS, cost += Norm(res);
    cost += std::log(1.0 + Norm(res) / lbfgs_dat->robust_nu);
  }

  // normalize cost by number of baselines
  return cost / double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline);
}

void ScalarGradient(double* unknowns, double* gradient, int n_unknowns,
                    void* extra_data) {
  assert(extra_data);
  const lbfgs_fulljones_data* lbfgs_dat =
      reinterpret_cast<lbfgs_fulljones_data*>(extra_data);
  const std::size_t n_sol = lbfgs_dat->n_antennas * lbfgs_dat->n_directions;
  std::vector<SolverBase::DComplex> solutions(n_sol);
  // Map unknowns into DComplex vector
  // Combine the real and imaginary values that are in the first and second half
  // of the unknowns array, respectively
  std::transform(unknowns, unknowns + n_sol, unknowns + n_sol,
                 solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });

  std::vector<SolverBase::DComplex> complex_gradient(
      n_sol, SolverBase::DComplex(0.0, 0.0));

  const double inv_baselines =
      1.0 / (double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline));
  for (size_t vis_index = lbfgs_dat->start_baseline;
       vis_index < lbfgs_dat->end_baseline; ++vis_index) {
    aocommon::MatrixComplexDouble2x2 res(
        lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
    // loop for residual calculation
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const SolverBase::DComplex* g_1 =
          &solutions[(antenna1 * lbfgs_dat->n_directions + d)];
      const SolverBase::DComplex* g_2 =
          &solutions[(antenna2 * lbfgs_dat->n_directions + d)];
      res -= model_data * (g_1[0] * std::conj(g_2[0]));
    }
    // Scale factor for robust gradient, divided by number_of_baselines
    // For LS cost, scale_factor=2.0/number_of_baselines
    const double scale_factor =
        -2.0 * inv_baselines /
        (lbfgs_dat->robust_nu + Norm(res));  //-ve for -ve grad direction
    // loop for grad calculation
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      // C
      const aocommon::MatrixComplexDouble2x2 model_data(
          lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const SolverBase::DComplex* g_1 =
          &solutions[(antenna1 * lbfgs_dat->n_directions + d)];
      const SolverBase::DComplex* g_2 =
          &solutions[(antenna2 * lbfgs_dat->n_directions + d)];

      // trace(C R^H)
      const SolverBase::DComplex traceCR = model_data[0] * std::conj(res[0]) +
                                           model_data[1] * std::conj(res[1]) +
                                           model_data[2] * std::conj(res[2]) +
                                           model_data[3] * std::conj(res[3]);

      // grad for antenna 1 = - g_2^H trace(C R^H)
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d)] +=
          scale_factor * std::conj(g_2[0]) * traceCR;
      // grad for antenna 2 = - g_1^H trace(C^H R)
      complex_gradient[(antenna2 * lbfgs_dat->n_directions + d)] +=
          scale_factor * std::conj(g_1[0] * traceCR);
    }
  }

  // Copy DComplex vector into gradient (real and imaginary parts separately)
  std::transform(complex_gradient.begin(), complex_gradient.end(), gradient,
                 [](const SolverBase::DComplex& z) { return z.real(); });
  // Note : return conjugate of gradient
  std::transform(complex_gradient.begin(), complex_gradient.end(),
                 gradient + n_sol,
                 [](const SolverBase::DComplex& z) { return -z.imag(); });
}

void PerformIterationFull(size_t ch_block,
                          const SolveData::ChannelBlockData& cb_data,
                          const std::vector<SolverBase::DComplex>& solutions,
                          SolutionSpan& next_solutions, std::size_t n_antennas,
                          std::size_t n_directions, double robust_nu,
                          std::size_t max_iter, std::size_t history_size,
                          std::size_t minibatches, persistent_data_t& pt) {
  lbfgs_fulljones_data t(cb_data, n_antennas, n_directions, robust_nu);
  const std::size_t n_solutions_2 = n_antennas * n_directions * 4;
  const std::size_t n_solutions = n_solutions_2 * 2;
  assert(n_solutions == static_cast<std::size_t>(pt.m));
  assert(cb_data.NVisibilities() == static_cast<std::size_t>(pt.nlen));

  xt::xtensor<double, 1> d_storage = LBFGSSolver::SplitSolutions(solutions);

  const std::size_t batch = (minibatches > 1 ? std::rand() % minibatches : 0);
  t.start_baseline = pt.offsets[batch];
  t.end_baseline = pt.offsets[batch] + pt.lengths[batch];
  // cannot work with just one baseline
  assert(t.end_baseline > t.start_baseline);
  lbfgs_fit(FullCost, FullGradient, d_storage.data(), n_solutions, max_iter,
            history_size, (void*)&t, &pt);

  LBFGSSolver::MergeSolutions(next_solutions, ch_block, d_storage);
}

void PerformIterationDiagonal(
    size_t ch_block, const SolveData::ChannelBlockData& cb_data,
    const std::vector<SolverBase::DComplex>& solutions,
    SolutionSpan& next_solutions, std::size_t n_antennas,
    std::size_t n_directions, double robust_nu, std::size_t max_iter,
    std::size_t history_size, std::size_t minibatches, persistent_data_t& pt) {
  lbfgs_fulljones_data t(cb_data, n_antennas, n_directions, robust_nu);
  const std::size_t n_solutions_2 = n_antennas * n_directions * 2;
  const std::size_t n_solutions = n_solutions_2 * 2;
  assert(n_solutions == static_cast<std::size_t>(pt.m));
  assert(cb_data.NVisibilities() == static_cast<std::size_t>(pt.nlen));

  xt::xtensor<double, 1> d_storage = LBFGSSolver::SplitSolutions(solutions);

  const std::size_t batch = (minibatches > 1 ? std::rand() % minibatches : 0);
  t.start_baseline = pt.offsets[batch];
  t.end_baseline = pt.offsets[batch] + pt.lengths[batch];
  // cannot work with just one baseline
  assert(t.end_baseline > t.start_baseline);
  lbfgs_fit(DiagonalCost, DiagonalGradient, d_storage.data(), n_solutions,
            max_iter, history_size, (void*)&t, &pt);

  LBFGSSolver::MergeSolutions(next_solutions, ch_block, d_storage);
}

void PerformIterationScalar(size_t ch_block,
                            const SolveData::ChannelBlockData& cb_data,
                            const std::vector<SolverBase::DComplex>& solutions,
                            SolutionSpan& next_solutions,
                            std::size_t n_antennas, std::size_t n_directions,
                            double robust_nu, std::size_t max_iter,
                            std::size_t history_size, std::size_t minibatches,
                            persistent_data_t& pt) {
  lbfgs_fulljones_data t(cb_data, n_antennas, n_directions, robust_nu);
  const std::size_t n_solutions_2 = n_antennas * n_directions;
  const std::size_t n_solutions = n_solutions_2 * 2;
  assert(n_solutions == static_cast<std::size_t>(pt.m));
  assert(cb_data.NVisibilities() == static_cast<std::size_t>(pt.nlen));

  xt::xtensor<double, 1> d_storage = LBFGSSolver::SplitSolutions(solutions);

  const std::size_t batch = (minibatches > 1 ? std::rand() % minibatches : 0);
  t.start_baseline = pt.offsets[batch];
  t.end_baseline = pt.offsets[batch] + pt.lengths[batch];
  // cannot work with just one baseline
  assert(t.end_baseline > t.start_baseline);
  lbfgs_fit(ScalarCost, ScalarGradient, d_storage.data(), n_solutions, max_iter,
            history_size, (void*)&t, &pt);

  LBFGSSolver::MergeSolutions(next_solutions, ch_block, d_storage);
}
}  // anonymous namespace

xt::xtensor<double, 1> LBFGSSolver::SplitSolutions(
    const std::vector<SolverBase::DComplex>& solutions) {
  const auto solutions_view = xt::adapt(solutions);
  const size_t n_solutions = solutions.size();

  xt::xtensor<double, 1> d_storage({n_solutions * 2},
                                   xt::layout_type::row_major);
  xt::view(d_storage, xt::range(0, n_solutions)) = xt::real(solutions_view);
  xt::view(d_storage, xt::range(n_solutions, n_solutions * 2)) =
      xt::imag(solutions_view);
  return d_storage;
}

void LBFGSSolver::MergeSolutions(SolutionSpan& next_solutions, size_t ch_block,
                                 const xt::xtensor<double, 1>& d_storage) {
  auto flat_next_solutions = xt::flatten(
      xt::view(next_solutions, ch_block, xt::all(), xt::all(), xt::all()));
  const size_t n_solutions = flat_next_solutions.size();

  assert(d_storage.size() == n_solutions * 2);
  auto d_storage_real = xt::view(d_storage, xt::range(0, n_solutions));
  auto d_storage_imag =
      xt::view(d_storage, xt::range(n_solutions, n_solutions * 2));

  flat_next_solutions = xt::vectorize([](double r, double i) {
    return SolverBase::DComplex(r, i);
  })(d_storage_real, d_storage_imag);
}

LBFGSSolver::SolveResult LBFGSSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  std::vector<persistent_data_t> persistent_data(NChannelBlocks());

  PrepareConstraints();

  SolveResult result;

  SolverMode mode(GetMode());

  SolutionTensor next_solutions_tensor(
      {NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()});
  SolutionSpan next_solutions = aocommon::xt::CreateSpan(next_solutions_tensor);

  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    const SolveData::ChannelBlockData& channel_block_data =
        data.ChannelBlock(ch_block);
    // initialize with
    // minibatches, size of parameters, size of data (baselines), history size,
    // 1 thread
    lbfgs_persist_init(
        &persistent_data[ch_block], GetMinibatches(),
        NAntennas() * NSolutions() * NSolutionPolarizations() * 2,
        channel_block_data.NVisibilities(), GetHistorySize(), 1);
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(GetMaxIterations());

  do {
    aocommon::ParallelFor<size_t> loop(GetNThreads());
    switch (mode) {
      case LBFGSSolver::kFull:
        MakeSolutionsFinite4Pol(solutions);
        loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t /*thread*/) {
          PerformIterationFull(ch_block, data.ChannelBlock(ch_block),
                               solutions[ch_block], next_solutions, NAntennas(),
                               NSolutions(), GetRobustDOF(), GetMaxIter(),
                               GetHistorySize(), GetMinibatches(),
                               persistent_data[ch_block]);
        });
        break;
      case LBFGSSolver::kDiagonal:
        MakeSolutionsFinite2Pol(solutions);
        loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t /*thread*/) {
          PerformIterationDiagonal(ch_block, data.ChannelBlock(ch_block),
                                   solutions[ch_block], next_solutions,
                                   NAntennas(), NSolutions(), GetRobustDOF(),
                                   GetMaxIter(), GetHistorySize(),
                                   GetMinibatches(), persistent_data[ch_block]);
        });
        break;
      case LBFGSSolver::kScalar:
        MakeSolutionsFinite1Pol(solutions);
        loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t /*thread*/) {
          PerformIterationScalar(ch_block, data.ChannelBlock(ch_block),
                                 solutions[ch_block], next_solutions,
                                 NAntennas(), NSolutions(), GetRobustDOF(),
                                 GetMaxIter(), GetHistorySize(),
                                 GetMinibatches(), persistent_data[ch_block]);
        });
        break;
      default:
        assert(false);
    }

    if (stat_stream) {
      (*stat_stream) << iteration << '\t';
    }

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes);

    if (stat_stream) {
      (*stat_stream) << step_magnitudes.back() << '\t' << avg_squared_diff
                     << '\n';
    }
    iteration++;

    has_previously_converged = has_converged || has_previously_converged;

  } while (!ReachedStoppingCriterion(iteration, has_converged,
                                     constraints_satisfied, step_magnitudes));

  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    lbfgs_persist_clear(&persistent_data[ch_block]);
  }

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged solves can be distinguished from converged ones.
  if (has_converged && constraints_satisfied)
    result.iterations = iteration;
  else
    result.iterations = iteration + 1;
  return result;
}
#endif /* HAVE_LIBDIRAC */
}  // namespace ddecal
}  // namespace dp3
