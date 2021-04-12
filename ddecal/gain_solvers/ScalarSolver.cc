// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ScalarSolver.h"

#include "SolverBuffer.h"

#include "../../base/DPBuffer.h"
#include "../linear_solvers/LLSSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <iomanip>
#include <iostream>
#include <boost/make_unique.hpp>

namespace dp3 {
namespace base {

ScalarSolver::SolveResult ScalarSolver::Solve(
    const SolverBuffer& solver_buffer,
    std::vector<std::vector<DComplex>>& solutions, double time,
    std::ostream* stat_stream) {
  assert(solutions.size() == NChannelBlocks());

  const size_t n_times = solver_buffer.NTimes();

  for (Constraint* c : GetConstraints()) c->PrepareIteration(false, 0, false);

  std::vector<std::vector<DComplex>> next_solutions(NChannelBlocks());

  SolveResult result;
  result.results.resize(GetConstraints().size());

  // Model matrix ant x [N x D] and visibility vector ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<Matrix>> g_times_cs(NChannelBlocks());
  std::vector<std::vector<Matrix>> vs(NChannelBlocks());
  for (size_t chBlock = 0; chBlock != NChannelBlocks(); ++chBlock) {
    next_solutions[chBlock].resize(NDirections() * NAntennas());
    const size_t channel_index_start = FirstChannel(chBlock);
    const size_t channel_index_end = FirstChannel(chBlock + 1);
    const size_t cur_channel_block_size =
        channel_index_end - channel_index_start;
    g_times_cs[chBlock].reserve(NAntennas());
    vs[chBlock].reserve(NAntennas());

    for (size_t ant = 0; ant != NAntennas(); ++ant) {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set
      // to 0.
      const size_t m = NAntennas() * n_times * cur_channel_block_size * 4;
      const size_t n = NDirections();
      const size_t nrhs = 1;
      g_times_cs[chBlock].emplace_back(m, n);
      vs[chBlock].emplace_back(std::max(m, n), nrhs);
    }
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  size_t constrained_iterations = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;
  bool has_stalled = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(GetMaxIterations());

  double avg_squared_diff = 1.0E4;

  do {
    MakeSolutionsFinite1Pol(solutions);

    ParallelFor<size_t> loop(GetNThreads());
    loop.Run(0, NChannelBlocks(), [&](size_t chBlock, size_t /*thread*/) {
      PerformIteration(solver_buffer, chBlock, g_times_cs[chBlock], vs[chBlock],
                       solutions[chBlock], next_solutions[chBlock],
                       double(iteration + 1) / GetMaxIterations(),
                       avg_squared_diff);
    });

    Step(solutions, next_solutions);

    if (stat_stream) {
      (*stat_stream) << iteration << '\t';
    }

    constraints_satisfied = true;

    for (size_t i = 0; i != GetConstraints().size(); ++i) {
      Constraint* c = GetConstraints()[i];
      // PrepareIteration() might change Satisfied(), and since we always want
      // to iterate at least once more when a constraint is not yet satisfied,
      // we evaluate Satisfied() before preparing.
      constraints_satisfied = c->Satisfied() && constraints_satisfied;
      c->PrepareIteration(has_previously_converged, iteration,
                          iteration + 1 >= GetMaxIterations());
      result.results[i] = c->Apply(next_solutions, time, stat_stream);
    }

    if (!constraints_satisfied) constrained_iterations = iteration + 1;

    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes, 1);

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
  if ((!has_converged || !constraints_satisfied) && !has_stalled)
    result.iterations = iteration + 1;
  else
    result.iterations = iteration;
  result.constraint_iterations = constrained_iterations;
  return result;
}

void ScalarSolver::PerformIteration(const SolverBuffer& solver_buffer,
                                    size_t channel_block_index,
                                    std::vector<Matrix>& g_times_cs,
                                    std::vector<Matrix>& vs,
                                    const std::vector<DComplex>& solutions,
                                    std::vector<DComplex>& next_solutions,
                                    double iterationfraction,
                                    double solverprecision) {
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    g_times_cs[ant].SetZero();
    vs[ant].SetZero();
  }

  const size_t channel_index_start = FirstChannel(channel_block_index);
  const size_t channel_index_end = FirstChannel(channel_block_index + 1);
  const size_t cur_channel_block_size = channel_index_end - channel_index_start;
  const size_t n_times = solver_buffer.NTimes();

  // The following loop fills the matrices for all antennas
  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    std::vector<const Complex*> model_ptrs(NDirections());
    for (size_t baseline = 0; baseline != NBaselines(); ++baseline) {
      size_t antenna1 = AntennaIndex1(baseline);
      size_t antenna2 = AntennaIndex2(baseline);
      if (antenna1 != antenna2) {
        Matrix& g_times_c1 = g_times_cs[antenna1];
        Matrix& v1 = vs[antenna1];
        Matrix& g_times_c2 = g_times_cs[antenna2];
        Matrix& v2 = vs[antenna2];
        for (size_t d = 0; d != NDirections(); ++d) {
          model_ptrs[d] = solver_buffer.ModelDataPointer(
              time_index, d, baseline, channel_index_start);
        }
        const Complex* data_ptr = solver_buffer.DataPointer(
            time_index, baseline, channel_index_start);
        const size_t p1_to_p2[4] = {0, 2, 1, 3};
        for (size_t ch = 0; ch < cur_channel_block_size; ++ch) {
          const size_t data_index1 =
              ch + (time_index + antenna1 * n_times) * cur_channel_block_size;
          const size_t data_index2 =
              ch + (time_index + antenna2 * n_times) * cur_channel_block_size;
          for (size_t p1 = 0; p1 != 4; ++p1) {
            size_t p2 = p1_to_p2[p1];
            for (size_t d = 0; d != NDirections(); ++d) {
              std::complex<double> predicted = *model_ptrs[d];

              size_t sol_index1 = antenna1 * NDirections() + d;
              size_t sol_index2 = antenna2 * NDirections() + d;
              g_times_c2(data_index1 * 4 + p1, d) = std::conj(
                  solutions[sol_index1] * predicted);  // using a* b* = (ab)*
              g_times_c1(data_index2 * 4 + p2, d) =
                  std::conj(solutions[sol_index2]) * predicted;

              ++model_ptrs[d];  // Goto the next polarization of this 2x2
                                // matrix.
            }
            v1(data_index2 * 4 + p2, 0) = *data_ptr;
            v2(data_index1 * 4 + p1, 0) = std::conj(*data_ptr);
            ++data_ptr;  // Goto the next polarization of this 2x2 matrix.
          }
        }
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  const size_t m = NAntennas() * n_times * cur_channel_block_size * 4;
  const size_t n = NDirections();
  const size_t nrhs = 1;
  std::unique_ptr<LLSSolver> solver =
      CreateLLSSolver(m, n, nrhs, iterationfraction, solverprecision);

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    // solve x^H in [g C] x^H  = v
    std::vector<Complex> x0(NDirections());
    for (size_t d = 0; d != NDirections(); ++d) {
      x0[d] = solutions[(ant * NDirections() + d)];
    }
    bool success =
        solver->Solve(g_times_cs[ant].data(), vs[ant].data(), x0.data());
    Matrix& x = vs[ant];
    if (success && x(0, 0) != Complex(0.0, 0.0)) {
      for (size_t d = 0; d != NDirections(); ++d)
        next_solutions[ant * NDirections() + d] = x(d, 0);
    } else {
      for (size_t d = 0; d != NDirections(); ++d)
        next_solutions[ant * NDirections() + d] =
            std::numeric_limits<double>::quiet_NaN();
    }
  }
}

}  // namespace base
}  // namespace dp3
