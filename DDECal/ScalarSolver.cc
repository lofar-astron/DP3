// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "QRSolver.h"
#include "ScalarSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <iomanip>
#include <iostream>

namespace DP3 {
namespace DPPP {

ScalarSolver::SolveResult ScalarSolver::Solve(
    const std::vector<Complex*>& unweighted_data,
    const std::vector<float*>& weights,
    std::vector<std::vector<Complex*> >&& unweighted_model_data,
    std::vector<std::vector<DComplex> >& solutions, double time,
    std::ostream* stat_stream) {
  const size_t n_times = unweighted_data.size();

  buffer_.AssignAndWeight(unweighted_data, weights,
                          std::move(unweighted_model_data));

  for (size_t i = 0; i != constraints_.size(); ++i)
    constraints_[i]->PrepareIteration(false, 0, false);

  std::vector<std::vector<DComplex> > next_solutions(n_channel_blocks_);

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
  std::vector<std::vector<Matrix> > g_times_cs(n_channel_blocks_);
  std::vector<std::vector<Matrix> > vs(n_channel_blocks_);
  for (size_t chBlock = 0; chBlock != n_channel_blocks_; ++chBlock) {
    next_solutions[chBlock].resize(n_directions_ * n_antennas_);
    const size_t channel_index_start =
                     chBlock * n_channels_ / n_channel_blocks_,
                 channel_index_end =
                     (chBlock + 1) * n_channels_ / n_channel_blocks_,
                 cur_channel_block_size =
                     channel_index_end - channel_index_start;
    g_times_cs[chBlock].resize(n_antennas_);
    vs[chBlock].resize(n_antennas_);

    for (size_t ant = 0; ant != n_antennas_; ++ant) {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set
      // to 0.
      size_t m = n_antennas_ * n_times * cur_channel_block_size * 4,
             n = n_directions_, nrhs = 1;
      g_times_cs[chBlock][ant] = Matrix(m, n);
      vs[chBlock][ant] = Matrix(std::max(m, n), nrhs);
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

  do {
    MakeSolutionsFinite1Pol(solutions);

    ParallelFor<size_t> loop(n_threads_);
    loop.Run(0, n_channel_blocks_, [&](size_t chBlock, size_t /*thread*/) {
      PerformIteration(chBlock, g_times_cs[chBlock], vs[chBlock],
                       solutions[chBlock], next_solutions[chBlock]);
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

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes, 1);
    if (stat_stream) {
      (*stat_stream) << step_magnitudes.back() << '\t' << avg_squared_diff
                     << '\n';
    }
    iteration++;

    has_previously_converged = has_converged || has_previously_converged;

    if (detect_stalling_ && constraints_satisfied)
      has_stalled = DetectStall(iteration, step_magnitudes);

  } while (iteration < max_iterations_ &&
           (!has_converged || !constraints_satisfied) && !has_stalled);

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged iterations can be distinguished from converged ones.
  if ((!has_converged || !constraints_satisfied) && !has_stalled)
    result.iterations = iteration + 1;
  else
    result.iterations = iteration;
  result.constraint_iterations = constrained_iterations;
  return result;
}

void ScalarSolver::PerformIteration(size_t channel_block_index,
                                    std::vector<Matrix>& g_times_cs,
                                    std::vector<Matrix>& vs,
                                    const std::vector<DComplex>& solutions,
                                    std::vector<DComplex>& next_solutions) {
  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    g_times_cs[ant].SetZero();
    vs[ant].SetZero();
  }

  const size_t channel_index_start =
      channel_block_index * n_channels_ / n_channel_blocks_;
  const size_t channel_index_end =
      (channel_block_index + 1) * n_channels_ / n_channel_blocks_;
  const size_t cur_channel_block_size = channel_index_end - channel_index_start,
               n_times = buffer_.Data().size();

  // The following loop fills the matrices for all antennas
  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    std::vector<const Complex*> model_ptrs(n_directions_);
    for (size_t baseline = 0; baseline != ant1_.size(); ++baseline) {
      size_t antenna1 = ant1_[baseline];
      size_t antenna2 = ant2_[baseline];
      if (antenna1 != antenna2) {
        Matrix& g_times_c1 = g_times_cs[antenna1];
        Matrix& v1 = vs[antenna1];
        Matrix& g_times_c2 = g_times_cs[antenna2];
        Matrix& v2 = vs[antenna2];
        for (size_t d = 0; d != n_directions_; ++d) {
          model_ptrs[d] =
              &buffer_.ModelData()[time_index][d][(channel_index_start +
                                                   baseline * n_channels_) *
                                                  4];
        }
        const Complex* data_ptr =
            &buffer_.Data()[time_index]
                           [(channel_index_start + baseline * n_channels_) * 4];
        const size_t p1_to_p2[4] = {0, 2, 1, 3};
        for (size_t ch = channel_index_start; ch != channel_index_end; ++ch) {
          const size_t data_index1 = ch - channel_index_start +
                                     (time_index + antenna1 * n_times) *
                                         cur_channel_block_size,
                       data_index2 = ch - channel_index_start +
                                     (time_index + antenna2 * n_times) *
                                         cur_channel_block_size;
          for (size_t p1 = 0; p1 != 4; ++p1) {
            size_t p2 = p1_to_p2[p1];
            for (size_t d = 0; d != n_directions_; ++d) {
              std::complex<double> predicted = *model_ptrs[d];

              size_t sol_index1 = antenna1 * n_directions_ + d;
              size_t sol_index2 = antenna2 * n_directions_ + d;
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
  size_t m = n_antennas_ * n_times * cur_channel_block_size * 4;
  size_t n = n_directions_, nrhs = 1;
  QRSolver solver(m, n, nrhs);
  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    // solve x^H in [g C] x^H  = v
    bool success = solver.Solve(g_times_cs[ant].data(), vs[ant].data());
    Matrix& x = vs[ant];
    if (success && x(0, 0) != Complex(0.0, 0.0)) {
      for (size_t d = 0; d != n_directions_; ++d)
        next_solutions[ant * n_directions_ + d] = x(d, 0);
    } else {
      for (size_t d = 0; d != n_directions_; ++d)
        next_solutions[ant * n_directions_ + d] =
            std::numeric_limits<double>::quiet_NaN();
    }
  }
}

}  // namespace DPPP
}  // namespace DP3
