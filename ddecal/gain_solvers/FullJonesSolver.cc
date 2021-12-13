// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "FullJonesSolver.h"

#include "../linear_solvers/QRSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

#include <boost/make_unique.hpp>

#include <iostream>

namespace dp3 {
namespace ddecal {

FullJonesSolver::SolveResult FullJonesSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  // This algorithm is basically the same as the scalar algorithm,
  // but visibility values are extended to 2x2 matrices and concatenated
  // in the matrix equations as block matrices. One difference is that
  // order of operations are important because of the non-commutativity of
  // matrix multiplication, as well as that A^H B^H = [BA]^H.
  //
  // The approach:
  // First we pre-apply the left-hand solutions to the model to make JM. Each
  // 2x2 coherence matrix Ji is matrix-multied by the lh solutions, for all
  // directions, and visibilities (times x channels).
  //   JMi = Ji Mi
  // These are stacked in matrix JM :
  //        JM0_d0 JM1_d0 ...
  //   JM = JM0_d1 JM1_d1
  //        ...
  // such that JM is a (2D) rows x (2N) col matrix, N=nvis, D=ndir.
  // The solved 2D x 2 solution matrix is similarly formed with the solution
  // values:
  //       ( J0 )
  //   J = ( J1 )
  //       ( .. )
  // And the 2N x 2 visibility matrix as well:
  //       ( V0 )
  //   V = ( V1 )
  //       ( .. )
  // And we solve the equation:
  //   'JM' J^H = V
  // With dimensions:
  //   [ 2N x 2D ] [ 2D x 2 ] = [ 2N x 2 ]

  assert(solutions.size() == NChannelBlocks());

  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(NChannelBlocks());

  SolveResult result;

  // Dimensions for each channelblock:
  // - Model matrix: n_antennas x [2N x 2D]
  // - Visibility matrix: n_antennas x [2N x 2]
  // The following loop allocates all structures:
  std::vector<std::vector<Matrix>> g_times_cs(NChannelBlocks());
  std::vector<std::vector<Matrix>> vs(NChannelBlocks());
  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    const SolveData::ChannelBlockData& channel_block_data =
        data.ChannelBlock(ch_block);
    next_solutions[ch_block].resize(NDirections() * NAntennas() * 4);
    g_times_cs[ch_block].reserve(NAntennas());
    vs[ch_block].reserve(NAntennas());

    for (size_t ant = 0; ant != NAntennas(); ++ant) {
      // Model matrix [2N x 2D] and visibility matrix [2N x 2]
      const size_t m = channel_block_data.NAntennaVisibilities(ant) * 2;
      const size_t n = NDirections() * 2;
      const size_t n_rhs = 2;
      g_times_cs[ch_block].emplace_back(m, n);
      vs[ch_block].emplace_back(std::max(m, n), n_rhs);
    }
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
    loop.Run(0, NChannelBlocks(),
             [&](size_t ch_block, [[maybe_unused]] size_t thread) {
               PerformIteration(data.ChannelBlock(ch_block),
                                g_times_cs[ch_block], vs[ch_block],
                                solutions[ch_block], next_solutions[ch_block]);
             });

    Step(solutions, next_solutions);

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

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged solves can be distinguished from converged ones.
  if (has_converged && constraints_satisfied)
    result.iterations = iteration;
  else
    result.iterations = iteration + 1;
  return result;
}

void FullJonesSolver::PerformIteration(
    const SolveData::ChannelBlockData& cb_data, std::vector<Matrix>& g_times_cs,
    std::vector<Matrix>& vs, const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  using aocommon::MC2x2;

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    g_times_cs[ant].SetZero();
    vs[ant].SetZero();
  }

  // Solving antenna2 uses visibilities of the 'normal' correlation
  // antenna1 x antenna2^H. The equation is:
  //   J_1 M J_2^H = V
  // Since in this equation antenna2 is solved, the solve matrices are
  // called g_times_c2 and v2. The index into these matrices is
  // ant2_positions[antenna2], which is increased in each iteration.
  //
  // Solving antenna1 uses visibilities of correlation antenna2 x antenna1^H.
  // It requires an extra Herm transpose on M and V. The equation is:
  //   J_2 M^H J_1^H = V^H
  // The relevant matrices and index are called g_times_c1, v1 and
  // ant1_positions[antenna1].

  // Using a single loop that fills both vs and g_times_cs is possible and it
  // would remove the duplicate code for ant1_positions and ant2_positions.
  // However, a performance test using the solvers/full_jones unit test
  // showed that using separate loops for each set of visibilities is slightly
  // faster. Using 1 thread (by adjusting kNThreads in SolverTester.h), the
  // total time spent in these loops was as follows:
  // One big loop: 27,3 seconds.
  // Separate loops: 24,6 seconds.

  // The following loop fills vs (for all antennas)
  std::vector<size_t> ant_positions(NAntennas(), 0);
  for (size_t vis_index = 0; vis_index != cb_data.NVisibilities();
       ++vis_index) {
    const aocommon::MC2x2F& data = cb_data.Visibility(vis_index);
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    size_t& a1pos = ant_positions[antenna1];
    size_t& a2pos = ant_positions[antenna2];

    Matrix& v1 = vs[antenna1];
    Matrix& v2 = vs[antenna2];
    for (size_t p = 0; p != 4; ++p) {
      v1(a1pos + (p % 2), p / 2) = std::conj(data[p]);
      // Perform the Hermitian transpose for v2 by flipping the p%2 and p/2
      // indices and omitting 'std::conj'.
      v2(a2pos + (p / 2), p % 2) = data[p];
    }

    a1pos += 2;
    a2pos += 2;
  }

  // The following loop fills g_times_cs (for all antennas)
  for (size_t d = 0; d != NDirections(); ++d) {
    ant_positions.assign(NAntennas(), 0);
    for (size_t vis_index = 0; vis_index != cb_data.NVisibilities();
         ++vis_index) {
      const size_t antenna1 = cb_data.Antenna1Index(vis_index);
      const size_t antenna2 = cb_data.Antenna2Index(vis_index);
      size_t& a1pos = ant_positions[antenna1];
      size_t& a2pos = ant_positions[antenna2];

      Matrix& g_times_c1 = g_times_cs[antenna1];
      Matrix& g_times_c2 = g_times_cs[antenna2];

      // Converting the solutions to std::complex<float> and using single
      // precision for the computation below, reduced the performance.
      // -> Keep using double precision until 'solutions' uses single precision.
      const MC2x2 model_data(cb_data.ModelVisibility(d, vis_index).Data());
      const MC2x2 solutions1(&solutions[(antenna1 * NDirections() + d) * 4]);
      const MC2x2 solutions2(&solutions[(antenna2 * NDirections() + d) * 4]);
      const MC2x2 g_times_c2_data = solutions1.Multiply(model_data);
      const MC2x2 g_times_c1_data = solutions2.MultiplyHerm(model_data);

      for (size_t p = 0; p != 4; ++p) {
        g_times_c2(a2pos + (p / 2), (d * 2) + (p % 2)) = g_times_c2_data[p];
        g_times_c1(a1pos + (p / 2), (d * 2) + (p % 2)) = g_times_c1_data[p];
      }

      a1pos += 2;
      a2pos += 2;
    }
  }

  // Compute the linear solution for each antenna.
  const size_t n = NDirections() * 2;
  const size_t n_rhs = 2;

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    const size_t m = cb_data.NAntennaVisibilities(ant) * 2;
    // TODO it would be nice to have a solver resize function to avoid too many
    // reallocations
    // TODO Create the solver using CreateLLSSolver, like the other solvers do.
    QRSolver solver(m, n, n_rhs);

    // solve x^H in [g C] x^H  = v
    bool success = solver.Solve(g_times_cs[ant].data(), vs[ant].data());
    Matrix& x = vs[ant];
    if (success && x(0, 0) != Complex(0.0, 0.0)) {
      for (size_t d = 0; d != NDirections(); ++d)
        for (size_t p = 0; p != 4; ++p) {
          // The conj transpose is also performed at this point (note swap of %
          // and /)
          next_solutions[(ant * NDirections() + d) * 4 + p] =
              std::conj(x(d * 2 + p % 2, p / 2));
        }
    } else {
      for (size_t i = 0; i != NDirections() * 4; ++i)
        next_solutions[ant * NDirections() * 4 + i] =
            std::numeric_limits<double>::quiet_NaN();
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
