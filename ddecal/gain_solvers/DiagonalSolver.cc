// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DiagonalSolver.h"
#include "SolveData.h"

#include "../linear_solvers/LLSSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <iomanip>
#include <iostream>

namespace dp3 {
namespace ddecal {

DiagonalSolver::SolveResult DiagonalSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  PrepareConstraints();
  SolveResult result;

  std::vector<std::vector<DComplex>> next_solutions(NChannelBlocks());
  for (std::vector<DComplex>& next_solution : next_solutions)
    next_solution.resize(NDirections() * NAntennas() * 2);

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

  // Using more threads wastes CPU and memory resources.
  const size_t n_threads = std::min(NChannelBlocks(), GetNThreads());
  // For each thread:
  // - Model matrix 2 x ant x [2N x D]
  // - Visibility vector 2 x ant x [2N]
  std::vector<std::vector<Matrix>> thread_g_times_cs(n_threads);
  std::vector<std::vector<std::vector<Complex>>> thread_vs(n_threads);

  aocommon::ParallelFor<size_t> loop(n_threads);
  do {
    MakeSolutionsFinite2Pol(solutions);

    loop.Run(0, NChannelBlocks(),
             [&](size_t ch_block, [[maybe_unused]] size_t thread) {
               const SolveData::ChannelBlockData& channel_block =
                   data.ChannelBlock(ch_block);

               std::vector<Matrix>& g_times_cs = thread_g_times_cs[thread];
               std::vector<std::vector<Complex>>& vs = thread_vs[thread];
               InitializeModelMatrix(channel_block, g_times_cs, vs);

               PerformIteration(channel_block, g_times_cs, vs,
                                solutions[ch_block], next_solutions[ch_block]);
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

void DiagonalSolver::PerformIteration(
    const SolveData::ChannelBlockData& cb_data, std::vector<Matrix>& g_times_cs,
    std::vector<std::vector<Complex>>& vs,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  for (size_t ant_and_pol = 0; ant_and_pol != NAntennas() * 2; ++ant_and_pol) {
    g_times_cs[ant_and_pol].SetZero();
    vs[ant_and_pol].assign(vs[ant_and_pol].size(), 0);
  }
  const size_t n_visibilities = cb_data.NVisibilities();

  // The following loop fills vs (for all antennas)
  std::vector<size_t> ant_positions(NAntennas() * 2, 0);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const aocommon::MC2x2F d = cb_data.Visibility(vis_index);
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    for (size_t p = 0; p != 4; ++p) {
      const size_t p1 = p / 2;
      const size_t p2 = p % 2;
      std::vector<Complex>& v1 = vs[antenna1 * 2 + p1];
      std::vector<Complex>& v2 = vs[antenna2 * 2 + p2];
      size_t& a1pos = ant_positions[antenna1 * 2 + p1];
      size_t& a2pos = ant_positions[antenna2 * 2 + p2];

      v1[a1pos] = d[p];
      v2[a2pos] = std::conj(d[p]);
      ++a1pos;
      ++a2pos;
    }
  }

  // The following loop fills g_times_cs (for all antennas)
  for (size_t d = 0; d != NDirections(); ++d) {
    ant_positions.assign(NAntennas() * 2, 0);

    for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
      size_t antenna1 = cb_data.Antenna1Index(vis_index);
      size_t antenna2 = cb_data.Antenna2Index(vis_index);
      const aocommon::MC2x2F predicted = cb_data.ModelVisibility(d, vis_index);

      for (size_t p = 0; p != 4; ++p) {
        const size_t p1 = p / 2;
        const size_t p2 = p % 2;

        Matrix& g_times_c1 = g_times_cs[antenna1 * 2 + p1];
        Matrix& g_times_c2 = g_times_cs[antenna2 * 2 + p2];
        size_t& a1pos = ant_positions[antenna1 * 2 + p1];
        size_t& a2pos = ant_positions[antenna2 * 2 + p2];
        const size_t sol_index1 = (antenna1 * NDirections() + d) * 2 + p1;
        const size_t sol_index2 = (antenna2 * NDirections() + d) * 2 + p2;

        g_times_c1(a1pos, d) =
            std::conj(Complex(solutions[sol_index2])) * predicted[p];
        g_times_c2(a2pos, d) = std::conj(Complex(solutions[sol_index1]) *
                                         predicted[p]);  // using a* b* = (ab)*
        ++a1pos;
        ++a2pos;
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  const size_t n = NDirections();
  const size_t nrhs = 1;

  std::vector<Complex> x0(NDirections());
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    for (size_t pol = 0; pol != 2; ++pol) {
      std::vector<Complex>& x = vs[ant * 2 + pol];
      const size_t m = x.size();
      // TODO it would be nice to have a solver resize function to avoid too
      // many reallocations
      std::unique_ptr<LLSSolver> solver = CreateLLSSolver(m, n, nrhs);
      // solve x^H in [g C] x^H  = v
      for (size_t d = 0; d != NDirections(); ++d) {
        x0[d] = solutions[(ant * NDirections() + d) * 2 + pol];
      }
      bool success =
          solver->Solve(g_times_cs[ant * 2 + pol].data(), x.data(), x0.data());
      if (success && x[0] != Complex(0.0, 0.0)) {
        for (size_t d = 0; d != NDirections(); ++d)
          next_solutions[(ant * NDirections() + d) * 2 + pol] = x[d];
      } else {
        for (size_t d = 0; d != NDirections(); ++d)
          next_solutions[(ant * NDirections() + d) * 2 + pol] =
              std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void DiagonalSolver::InitializeModelMatrix(
    const SolveData::ChannelBlockData& channel_block_data,
    std::vector<Matrix>& g_times_cs,
    std::vector<std::vector<Complex>>& vs) const {
  assert(g_times_cs.empty() == vs.empty());
  if (g_times_cs.empty()) {
    // Executed the first iteration only.
    g_times_cs.reserve(NAntennas() * 2);
    vs.reserve(NAntennas() * 2);
  } else {
    // Note clear() does not modify the capacity.
    // See https://en.cppreference.com/w/cpp/container/vector/clear
    g_times_cs.clear();
    vs.clear();
  }

  // Update the size of the model matrix and initialize them to zero.
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    // Model matrix [(2N) x D] and visibility vector [2N]
    const size_t n_visibilities = channel_block_data.NAntennaVisibilities(ant);
    const size_t m = n_visibilities * 2;
    const size_t n = NDirections();
    for (size_t p = 0; p != 2; ++p) {
      g_times_cs.emplace_back(m, n);
      vs.emplace_back(std::max(m, n));
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
