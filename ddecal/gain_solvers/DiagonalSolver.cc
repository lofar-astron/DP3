// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DiagonalSolver.h"

#include "SolverBuffer.h"

#include "../../base/DPBuffer.h"
#include "../linear_solvers/LLSSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <iomanip>
#include <iostream>

namespace dp3 {
namespace base {

DiagonalSolver::SolveResult DiagonalSolver::Solve(
    const SolverBuffer& solver_buffer,
    std::vector<std::vector<DComplex>>& solutions, double time,
    std::ostream* stat_stream) {
  const size_t n_times = solver_buffer.NTimes();

  for (size_t i = 0; i != constraints_.size(); ++i)
    constraints_[i]->PrepareIteration(false, 0, false);

  std::vector<std::vector<DComplex>> next_solutions(n_channel_blocks_);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != n_channel_blocks_) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  result.results.resize(constraints_.size());

  // Model matrix ant x [N x D] and visibility vector ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<Matrix>> g_times_cs(n_channel_blocks_);
  std::vector<std::vector<std::vector<Complex>>> vs(n_channel_blocks_);
  for (size_t chBlock = 0; chBlock != n_channel_blocks_; ++chBlock) {
    next_solutions[chBlock].resize(n_directions_ * n_antennas_ * 2);
    const size_t channel_index_start =
                     chBlock * n_channels_ / n_channel_blocks_,
                 channel_index_end =
                     (chBlock + 1) * n_channels_ / n_channel_blocks_,
                 cur_channel_block_size =
                     channel_index_end - channel_index_start;
    g_times_cs[chBlock].resize(n_antennas_ * 2);
    vs[chBlock].resize(n_antennas_ * 2);

    for (size_t ant = 0; ant != n_antennas_ * 2; ++ant) {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set
      // to 0.
      // X and Y polarizations are treated as two different antennas.
      size_t m = (n_antennas_ * 2) * n_times * cur_channel_block_size;
      size_t n = n_directions_;
      g_times_cs[chBlock][ant] = Matrix(m, n);
      vs[chBlock][ant].resize(std::max(m, n));
    }
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0, constrained_iterations = 0;
  bool has_converged = false, has_previously_converged = false,
       constraints_satisfied = false, has_stalled = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(max_iterations_);

  double avg_squared_diff = 1.0E4;

  do {
    MakeSolutionsFinite2Pol(solutions);

    ParallelFor<size_t> loop(n_threads_);
    loop.Run(0, n_channel_blocks_, [&](size_t chBlock, size_t /*thread*/) {
      PerformIteration(solver_buffer, chBlock, g_times_cs[chBlock], vs[chBlock],
                       solutions[chBlock], next_solutions[chBlock],
                       (double)(iteration + 1) / max_iterations_,
                       avg_squared_diff);
    });

    Step(solutions, next_solutions);

    if (stat_stream) {
      (*stat_stream) << iteration << '\t';
    }

    constraints_satisfied = true;

    for (size_t i = 0; i != constraints_.size(); ++i) {
      // PrepareIteration() might change Satisfied(), and since we always want
      // to iterate at least once more when a constraint is not yet satisfied,
      // we evaluate Satisfied() before preparing.
      constraints_satisfied =
          constraints_[i]->Satisfied() && constraints_satisfied;
      constraints_[i]->PrepareIteration(has_previously_converged, iteration,
                                        iteration + 1 >= max_iterations_);
      result.results[i] =
          constraints_[i]->Apply(next_solutions, time, stat_stream);
    }

    if (!constraints_satisfied) constrained_iterations = iteration + 1;

    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes, 2);
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

void DiagonalSolver::PerformIteration(const SolverBuffer& solver_buffer,
                                      size_t channel_block_index,
                                      std::vector<Matrix>& g_times_cs,
                                      std::vector<std::vector<Complex>>& vs,
                                      const std::vector<DComplex>& solutions,
                                      std::vector<DComplex>& next_solutions,
                                      double iterationfraction,
                                      double solverprecision) {
  for (size_t ant = 0; ant != n_antennas_ * 2; ++ant) {
    g_times_cs[ant].SetZero();
    std::fill(vs[ant].begin(), vs[ant].end(), 0.0);
  }

  const size_t channel_index_start =
      channel_block_index * n_channels_ / n_channel_blocks_;
  const size_t channel_index_end =
      (channel_block_index + 1) * n_channels_ / n_channel_blocks_;
  const size_t cur_channel_block_size = channel_index_end - channel_index_start;
  const size_t n_times = solver_buffer.NTimes();

  // The following loop fills the matrices for all antennas
  std::vector<const Complex*> model_ptrs(n_directions_);
  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    for (size_t baseline = 0; baseline != ant1_.size(); ++baseline) {
      size_t antenna1 = ant1_[baseline];
      size_t antenna2 = ant2_[baseline];
      if (antenna1 != antenna2) {
        for (size_t d = 0; d != n_directions_; ++d) {
          model_ptrs[d] = solver_buffer.ModelDataPointer(
              time_index, d, baseline, channel_index_start);
        }
        const Complex* data_ptr = solver_buffer.DataPointer(
            time_index, baseline, channel_index_start);
        for (size_t ch = channel_index_start; ch != channel_index_end; ++ch) {
          for (size_t p = 0; p != 4; ++p) {
            size_t p1 = p / 2;
            size_t p2 = p % 2;
            const size_t data_index1 =
                             ch - channel_index_start +
                             (time_index + (antenna1 * 2 + p1) * n_times) *
                                 cur_channel_block_size,
                         data_index2 =
                             ch - channel_index_start +
                             (time_index + (antenna2 * 2 + p2) * n_times) *
                                 cur_channel_block_size;
            Matrix& g_times_c1 = g_times_cs[antenna1 * 2 + p1];
            std::vector<Complex>& v1 = vs[antenna1 * 2 + p1];
            Matrix& g_times_c2 = g_times_cs[antenna2 * 2 + p2];
            std::vector<Complex>& v2 = vs[antenna2 * 2 + p2];

            for (size_t d = 0; d != n_directions_; ++d) {
              std::complex<double> predicted = *model_ptrs[d];

              size_t sol_index1 = (antenna1 * n_directions_ + d) * 2 + p1;
              size_t sol_index2 = (antenna2 * n_directions_ + d) * 2 + p2;
              g_times_c2(data_index1, d) = std::conj(
                  solutions[sol_index1] * predicted);  // using a* b* = (ab)*
              g_times_c1(data_index2, d) =
                  std::conj(solutions[sol_index2]) * predicted;

              ++model_ptrs[d];  // Goto the next polarization of this 2x2
                                // matrix.
            }
            v1[data_index2] = *data_ptr;
            v2[data_index1] = std::conj(*data_ptr);
            ++data_ptr;  // Goto the next polarization of this 2x2 matrix.
          }
        }
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  const size_t m = n_antennas_ * 2 * n_times * cur_channel_block_size;
  const size_t n = n_directions_;
  const size_t nrhs = 1;
  std::unique_ptr<LLSSolver> solver =
      LLSSolver::Make(lls_solver_type_, m, n, nrhs);
  solver->SetTolerance(
      calculateLLSTolerance(iterationfraction, solverprecision));
  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    for (size_t pol = 0; pol != 2; ++pol) {
      // solve x^H in [g C] x^H  = v
      std::vector<Complex> x0(n_directions_);
      for (size_t d = 0; d != n_directions_; ++d) {
        x0[d] = solutions[(ant * n_directions_ + d) * 2 + pol];
      }

      bool success = solver->Solve(g_times_cs[ant * 2 + pol].data(),
                                   vs[ant * 2 + pol].data(), x0.data());
      std::vector<Complex>& x = vs[ant * 2 + pol];
      if (success && x[0] != Complex(0.0, 0.0)) {
        for (size_t d = 0; d != n_directions_; ++d)
          next_solutions[(ant * n_directions_ + d) * 2 + pol] = x[d];
      } else {
        for (size_t d = 0; d != n_directions_; ++d)
          next_solutions[(ant * n_directions_ + d) * 2 + pol] =
              std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

}  // namespace base
}  // namespace dp3
