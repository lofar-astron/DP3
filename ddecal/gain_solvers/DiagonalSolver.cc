// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DiagonalSolver.h"

#include <algorithm>
#include <iomanip>
#include <ostream>

#include <aocommon/matrix2x2.h>
#include <aocommon/recursivefor.h>

#include <xtensor/xview.hpp>

#include "../linear_solvers/LLSSolver.h"

#include "SolveData.h"

namespace dp3 {
namespace ddecal {

DiagonalSolver::SolveResult DiagonalSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  PrepareConstraints();
  SolveResult result;

  SolutionTensor next_solutions(
      {NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()});

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
  // Max 8 threads are used in the outer loop, to limit the required memory.
  // TODO this could be smarter; e.g. check required vs available memory.
  const size_t n_outer_threads = std::min<size_t>(8, n_threads);

  // For each thread:
  // - Model matrix 2 x ant x [2N x D]
  // - Visibility vector 2 x ant x [2N]
  std::vector<std::vector<Matrix>> thread_g_times_cs(n_outer_threads);
  std::vector<std::vector<std::vector<Complex>>> thread_vs(n_outer_threads);

  // A RecursiveFor is started to allow nested parallelization over antennas
  // in PerformIteration.
  aocommon::RecursiveFor recursive_for;
  do {
    MakeSolutionsFinite2Pol(solutions);

    recursive_for.ConstrainedRun(
        0, NChannelBlocks(), n_outer_threads,
        [&](size_t start_block, size_t end_block, size_t thread_index) {
          for (size_t ch_block = start_block; ch_block != end_block;
               ++ch_block) {
            const SolveData::ChannelBlockData& channel_block =
                data.ChannelBlock(ch_block);

            std::vector<Matrix>& g_times_cs = thread_g_times_cs[thread_index];
            std::vector<std::vector<Complex>>& vs = thread_vs[thread_index];
            InitializeModelMatrix(channel_block, g_times_cs, vs);

            PerformIteration(ch_block, channel_block, g_times_cs, vs,
                             solutions[ch_block], next_solutions);
          }
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
    size_t ch_block, const SolveData::ChannelBlockData& cb_data,
    std::vector<Matrix>& g_times_cs, std::vector<std::vector<Complex>>& vs,
    const std::vector<DComplex>& solutions, SolutionTensor& next_solutions) {
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
  for (size_t s = 0; s != NSolutions(); ++s) {
    ant_positions.assign(NAntennas() * 2, 0);

    for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
      size_t antenna1 = cb_data.Antenna1Index(vis_index);
      size_t antenna2 = cb_data.Antenna2Index(vis_index);
      const aocommon::MC2x2F predicted = cb_data.ModelVisibility(s, vis_index);

      for (size_t p = 0; p != 4; ++p) {
        const size_t p1 = p / 2;
        const size_t p2 = p % 2;

        Matrix& g_times_c1 = g_times_cs[antenna1 * 2 + p1];
        Matrix& g_times_c2 = g_times_cs[antenna2 * 2 + p2];
        size_t& a1pos = ant_positions[antenna1 * 2 + p1];
        size_t& a2pos = ant_positions[antenna2 * 2 + p2];
        const size_t sol_index1 = (antenna1 * NSolutions() + s) * 2 + p1;
        const size_t sol_index2 = (antenna2 * NSolutions() + s) * 2 + p2;

        g_times_c1(a1pos, s) =
            std::conj(Complex(solutions[sol_index2])) * predicted[p];
        g_times_c2(a2pos, s) = std::conj(Complex(solutions[sol_index1]) *
                                         predicted[p]);  // using a* b* = (ab)*
        ++a1pos;
        ++a2pos;
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  const size_t n = NSolutions();
  const size_t nrhs = 1;

  aocommon::RecursiveFor::NestedRun(0, NAntennas() * 2, [&](size_t i) {
    const size_t ant = i / 2;
    const size_t pol = i % 2;
    const size_t m = cb_data.NAntennaVisibilities(ant) * 2;
    std::unique_ptr<LLSSolver> solver = CreateLLSSolver(m, n, nrhs);
    // solve x^H in [g C] x^H  = v
    std::vector<Complex> x0(NSolutions());
    for (size_t s = 0; s != NSolutions(); ++s) {
      x0[s] = solutions[(ant * NSolutions() + s) * 2 + pol];
    }
    std::vector<Complex>& x = vs[ant * 2 + pol];
    bool success =
        solver->Solve(g_times_cs[ant * 2 + pol].data(), x.data(), x0.data());
    if (success && x[0] != Complex(0.0, 0.0)) {
      for (size_t s = 0; s != NSolutions(); ++s)
        next_solutions(ch_block, ant, s, pol) = x[s];
    } else {
      xt::view(next_solutions, ch_block, ant, xt::all(), pol)
          .fill(std::numeric_limits<double>::quiet_NaN());
    }
  });
}

// Based on SolverBase::Matrix::Reset.
static void Reset(std::vector<SolverBase::Complex>& vector, size_t size) {
  // Minimize the number of elements to modify.
  if (size < vector.size()) {
    vector.resize(size);
    std::fill(vector.begin(), vector.end(), SolverBase::Complex(0.0, 0.0));
  } else {
    std::fill(vector.begin(), vector.end(), SolverBase::Complex(0.0, 0.0));
    vector.resize(size, SolverBase::Complex(0.0, 0.0));
  }
}

void DiagonalSolver::InitializeModelMatrix(
    const SolveData::ChannelBlockData& channel_block_data,
    std::vector<Matrix>& g_times_cs,
    std::vector<std::vector<Complex>>& vs) const {
  assert(g_times_cs.empty() == vs.empty());
  if (g_times_cs.empty()) {
    // Executed the first iteration only.
    g_times_cs.resize(NAntennas() * 2);
    vs.resize(NAntennas() * 2);
  }

  // Update the size of the model matrix and initialize them to zero.
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    // Model matrix [(2N) x D] and visibility vector [2N]
    const size_t n_visibilities = channel_block_data.NAntennaVisibilities(ant);
    const size_t m = n_visibilities * 2;
    const size_t n = NSolutions();
    for (size_t p = 0; p != 2; ++p) {
      g_times_cs[2 * ant + p].Reset(m, n);
      Reset(vs[2 * ant + p], std::max(m, n));
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
