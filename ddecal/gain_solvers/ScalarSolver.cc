// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ScalarSolver.h"
#include "SolveData.h"

#include "../linear_solvers/LLSSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/dynamicfor.h>
#include <xtensor/xview.hpp>

#include <iostream>

namespace dp3 {
namespace ddecal {

ScalarSolver::SolveResult ScalarSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  PrepareConstraints();
  SolveResult result;

  SolutionTensor next_solutions_tensor(
      {NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()});
  SolutionSpan next_solutions = aocommon::xt::CreateSpan(next_solutions_tensor);

  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(GetMaxIterations());

  double avg_squared_diff = 1.0E4;

  const size_t n_threads = aocommon::ThreadPool::GetInstance().NThreads();
  // For each thread:
  // - Model matrix ant x [N x D] and
  // - Visibility vector ant x [N x 1]
  std::vector<std::vector<Matrix>> thread_g_times_cs(n_threads);
  std::vector<std::vector<Matrix>> thread_vs(n_threads);

  aocommon::DynamicFor<size_t> loop;
  do {
    MakeSolutionsFinite1Pol(solutions);

    loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t thread) {
      const SolveData::ChannelBlockData& channel_block =
          data.ChannelBlock(ch_block);

      std::vector<Matrix>& g_times_cs = thread_g_times_cs[thread];
      std::vector<Matrix>& vs = thread_vs[thread];
      InitializeModelMatrix(channel_block, g_times_cs, vs);

      PerformIteration(ch_block, channel_block, g_times_cs, vs,
                       solutions[ch_block], next_solutions);
    });

    Step(solutions, next_solutions);

    if (stat_stream) {
      (*stat_stream) << iteration << '\t';
    }

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

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

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged iterations can be distinguished from converged ones.
  if (has_converged && constraints_satisfied)
    result.iterations = iteration;
  else
    result.iterations = iteration + 1;
  return result;
}

void ScalarSolver::PerformIteration(size_t ch_block,
                                    const SolveData::ChannelBlockData& cb_data,
                                    std::vector<Matrix>& g_times_cs,
                                    std::vector<Matrix>& vs,
                                    const std::vector<DComplex>& solutions,
                                    SolutionSpan& next_solutions) {
  const size_t n_visibilities = cb_data.NVisibilities();
  const size_t p1_to_p2[4] = {0, 2, 1, 3};

  // The following loop fills vs (for all antennas)
  std::vector<size_t> ant_positions(NAntennas(), 0);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    size_t antenna1 = cb_data.Antenna1Index(vis_index);
    size_t antenna2 = cb_data.Antenna2Index(vis_index);
    Matrix& v1 = vs[antenna1];
    Matrix& v2 = vs[antenna2];

    const aocommon::MC2x2F d = cb_data.Visibility(vis_index);
    size_t& a1pos = ant_positions[antenna1];
    size_t& a2pos = ant_positions[antenna2];
    for (size_t p1 = 0; p1 != 4; ++p1) {
      v1(a1pos * 4 + p1_to_p2[p1], 0) = d[p1];
      v2(a2pos * 4 + p1, 0) = std::conj(d[p1]);
    }
    ++a1pos;
    ++a2pos;
  }

  // The following loop fills g_times_cs (for all antennas)
  for (size_t s = 0; s != NSolutions(); ++s) {
    ant_positions.assign(NAntennas(), 0);
    for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
      size_t antenna1 = cb_data.Antenna1Index(vis_index);
      size_t antenna2 = cb_data.Antenna2Index(vis_index);
      Matrix& g_times_c1 = g_times_cs[antenna1];
      Matrix& g_times_c2 = g_times_cs[antenna2];

      size_t& a1pos = ant_positions[antenna1];
      size_t& a2pos = ant_positions[antenna2];
      for (size_t p1 = 0; p1 != 4; ++p1) {
        const size_t p2 = p1_to_p2[p1];
        const aocommon::MC2x2F predicted =
            cb_data.ModelVisibility(s, vis_index);

        const size_t sol_index1 = antenna1 * NSolutions() + s;
        const size_t sol_index2 = antenna2 * NSolutions() + s;
        g_times_c1(a1pos * 4 + p2, s) =
            std::conj(Complex(solutions[sol_index2])) * predicted[p1];
        g_times_c2(a2pos * 4 + p1, s) =
            std::conj(Complex(solutions[sol_index1]) *
                      predicted[p1]);  // using a* b* = (ab)*
      }
      ++a1pos;
      ++a2pos;
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  const size_t n = NSolutions();
  const size_t nrhs = 1;
  std::vector<Complex> x0(NSolutions());

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    const size_t m = cb_data.NAntennaVisibilities(ant) * 4;
    // TODO it would be nice to have a solver resize function to avoid too many
    // reallocations
    std::unique_ptr<LLSSolver> solver = CreateLLSSolver(m, n, nrhs);
    // solve x^H in [g C] x^H  = v
    for (size_t s = 0; s != NSolutions(); ++s) {
      x0[s] = solutions[(ant * NSolutions() + s)];
    }
    Matrix& x = vs[ant];
    bool success = solver->Solve(g_times_cs[ant].data(), x.data(), x0.data());
    if (success && x(0, 0) != Complex(0.0, 0.0)) {
      for (size_t s = 0; s != NSolutions(); ++s)
        next_solutions(ch_block, ant, s, 0) = x(s, 0);
    } else {
      xt::view(next_solutions, ch_block, ant, xt::all(), 0)
          .fill(std::numeric_limits<double>::quiet_NaN());
    }
  }
}

void ScalarSolver::InitializeModelMatrix(
    const SolveData::ChannelBlockData& channel_block_data,
    std::vector<Matrix>& g_times_cs, std::vector<Matrix>& vs) const {
  assert(g_times_cs.empty() == vs.empty());
  if (g_times_cs.empty()) {
    // Executed the first iteration only.
    g_times_cs.resize(NAntennas());
    vs.resize(NAntennas());
  }

  // Update the size of the model matrix and initialize them to zero.
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    // Model matrix [N x D] and visibility vector [N x 1]
    const size_t n_visibilities = channel_block_data.NAntennaVisibilities(ant);
    const size_t m = n_visibilities * 4;
    const size_t n = NSolutions();
    const size_t n_rhs = 1;

    g_times_cs[ant].Reset(m, n);
    vs[ant].Reset(std::max(m, n), n_rhs);
  }
}

}  // namespace ddecal
}  // namespace dp3
