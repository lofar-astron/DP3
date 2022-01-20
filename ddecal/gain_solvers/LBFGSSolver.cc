// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "LBFGSSolver.h"

#ifdef HAVE_LIBDIRAC
#include <Dirac.h>
#endif /* HAVE_LIBDIRAC */

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

#include <boost/make_unique.hpp>

#include <iostream>
#include <algorithm>

namespace dp3 {
namespace ddecal {

#ifdef HAVE_LIBDIRAC
namespace {
using aocommon::MC2x2;
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

double fulljones_cost(double* unknowns, int n_uknowns, void* extra_data) {
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
    MC2x2 res(lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const MC2x2 model_data(lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const MC2x2 J1(&solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const MC2x2 J2(&solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);
      const MC2x2 J1_M = J1.Multiply(model_data);
      const MC2x2 J1_M_J2H = J1_M.MultiplyHerm(J2);
      res -= J1_M_J2H;
    }
    // For LS, cost += Norm(res);
    cost += std::log(1.0 + Norm(res) / lbfgs_dat->robust_nu);
  }

  // normalize cost by number of baselines
  return cost / double(lbfgs_dat->end_baseline - lbfgs_dat->start_baseline);
}

void fulljones_gradient(double* unknowns, double* gradient, int n_unknowns,
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
    MC2x2 res(lbfgs_dat->cb_data.Visibility(vis_index));
    const size_t antenna1 = lbfgs_dat->cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = lbfgs_dat->cb_data.Antenna2Index(vis_index);
    // loop for residual calculation
#pragma GCC ivdep
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const MC2x2 model_data(lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const MC2x2 J1(&solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const MC2x2 J2(&solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);
      const MC2x2 J1_M = J1.Multiply(model_data);
      const MC2x2 J1_M_J2H = J1_M.MultiplyHerm(J2);
      res -= J1_M_J2H;
    }
    // Scale factor for robust gradient, divided by number_of_baselines
    // For LS cost, scale_factor=2.0/number_of_baselines
    const double scale_factor =
        -2.0 * inv_baselines /
        (lbfgs_dat->robust_nu + Norm(res));  //-ve for -ve grad direction
    // loop for grad calculation
    for (size_t d = 0; d != lbfgs_dat->n_directions; ++d) {
      const MC2x2 model_data(lbfgs_dat->cb_data.ModelVisibility(d, vis_index));
      const MC2x2 J1(&solutions[(antenna1 * lbfgs_dat->n_directions + d) * 4]);
      const MC2x2 J2(&solutions[(antenna2 * lbfgs_dat->n_directions + d) * 4]);

      const MC2x2 res_J2 = res.Multiply(J2);
      const MC2x2 res_J2_modH = res_J2.MultiplyHerm(model_data);
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4] +=
          scale_factor * res_J2_modH[0];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 1] +=
          scale_factor * res_J2_modH[1];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 2] +=
          scale_factor * res_J2_modH[2];
      complex_gradient[(antenna1 * lbfgs_dat->n_directions + d) * 4 + 3] +=
          scale_factor * res_J2_modH[3];

      const MC2x2 resH_J1 = res.HermThenMultiply(J1);
      const MC2x2 resH_J1_mod = resH_J1.Multiply(model_data);
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

void PerformIteration(const SolveData::ChannelBlockData& cb_data,
                      const std::vector<SolverBase::DComplex>& solutions,
                      std::vector<SolverBase::DComplex>& next_solutions,
                      std::size_t n_antennas, std::size_t n_directions,
                      double robust_nu, std::size_t max_iter,
                      std::size_t history_size, std::size_t minibatches,
                      persistent_data_t& pt) {
  lbfgs_fulljones_data t(cb_data, n_antennas, n_directions, robust_nu);
  const std::size_t n_data = cb_data.NVisibilities() * 8;
  const std::size_t n_solutions_2 = n_antennas * n_directions * 4;
  const std::size_t n_solutions = n_solutions_2 * 2;
  assert(n_solutions == static_cast<std::size_t>(pt.m));
  assert(n_data == static_cast<std::size_t>(pt.nlen) * 8);

  std::vector<double> d_storage(n_solutions);
  std::transform(solutions.begin(), solutions.end(), d_storage.begin(),
                 [](const SolverBase::DComplex& z) { return z.real(); });
  std::transform(solutions.begin(), solutions.end(),
                 d_storage.begin() + n_solutions_2,
                 [](const SolverBase::DComplex& z) { return z.imag(); });

  const std::size_t batch = (minibatches > 1 ? std::rand() % minibatches : 0);
  t.start_baseline = pt.offsets[batch];
  t.end_baseline = pt.offsets[batch] + pt.lengths[batch];
  // cannot work with just one baseline
  assert(t.end_baseline > t.start_baseline);
  lbfgs_fit(fulljones_cost, fulljones_gradient, d_storage.data(), n_solutions,
            max_iter, history_size, (void*)&t, &pt);

  // Map double vector to solutions
  std::transform(d_storage.begin(), d_storage.begin() + n_solutions_2,
                 d_storage.begin() + n_solutions_2, next_solutions.begin(),
                 [](double r, double i) { return SolverBase::DComplex(r, i); });
}
}  // anonymous namespace

LBFGSSolver::SolveResult LBFGSSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  std::vector<persistent_data_t> persistent_data(NChannelBlocks());

  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(NChannelBlocks());

  SolveResult result;

  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    const SolveData::ChannelBlockData& channel_block_data =
        data.ChannelBlock(ch_block);
    next_solutions[ch_block].resize(NDirections() * NAntennas() * 4);
    // initialize with
    // minibatches, size of parameters, size of data (baselines), history size,
    // 1 thread
    lbfgs_persist_init(&persistent_data[ch_block], GetMinibatches(),
                       NDirections() * NAntennas() * 8,
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
    MakeSolutionsFinite4Pol(solutions);

    aocommon::ParallelFor<size_t> loop(GetNThreads());
    loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t /*thread*/) {
      PerformIteration(data.ChannelBlock(ch_block), solutions[ch_block],
                       next_solutions[ch_block], NAntennas(), NDirections(),
                       GetRobustDOF(), GetMaxIter(), GetHistorySize(),
                       GetMinibatches(), persistent_data[ch_block]);
    });
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
