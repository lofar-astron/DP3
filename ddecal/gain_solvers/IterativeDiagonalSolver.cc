// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeDiagonalSolver.h"

#include "SolverBuffer.h"

#include "../../base/DPBuffer.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace dp3 {
namespace ddecal {

IterativeDiagonalSolver::SolveResult IterativeDiagonalSolver::Solve(
    const SolverBuffer& solver_buffer,
    std::vector<std::vector<DComplex>>& solutions, double time,
    std::ostream* stat_stream) {
  const size_t n_times = solver_buffer.NTimes();

  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(n_channel_blocks_);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != n_channel_blocks_) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  // Visibility vector chblock x time x [bl * chan * ncor]
  constexpr size_t n_cor = 4;
  constexpr size_t n_solution_pols = 2;
  std::vector<std::vector<std::vector<Complex>>> v_residual(n_channel_blocks_);
  // The following loop allocates all structures
  for (size_t chBlock = 0; chBlock != n_channel_blocks_; ++chBlock) {
    next_solutions[chBlock].resize(NDirections() * NAntennas() *
                                   n_solution_pols);
    const size_t channel_index_start =
                     chBlock * n_channels_ / n_channel_blocks_,
                 channel_index_end =
                     (chBlock + 1) * n_channels_ / n_channel_blocks_,
                 cur_channel_block_size =
                     channel_index_end - channel_index_start;
    v_residual[chBlock].resize(n_times);

    for (std::vector<Complex>& data : v_residual[chBlock]) {
      data.resize(ant1_.size() * cur_channel_block_size * n_cor);
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
  step_magnitudes.reserve(max_iterations_);

  do {
    MakeSolutionsFinite2Pol(solutions);

    ParallelFor<size_t> loop(n_threads_);
    loop.Run(0, n_channel_blocks_, [&](size_t chBlock, size_t /*thread*/) {
      PerformIteration(solver_buffer, chBlock, v_residual[chBlock],
                       solutions[chBlock], next_solutions[chBlock]);
    });

    Step(solutions, next_solutions);

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes, n_solution_pols);
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

void IterativeDiagonalSolver::PerformIteration(
    const SolverBuffer& solver_buffer, size_t channel_block_index,
    std::vector<std::vector<Complex>>& v_residual,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  const size_t channel_index_start =
      channel_block_index * n_channels_ / n_channel_blocks_;
  const size_t channel_index_end =
      (channel_block_index + 1) * n_channels_ / n_channel_blocks_;
  const size_t n_times = solver_buffer.NTimes();

  // Fill v_residual
  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    solver_buffer.CopyDataChannels(time_index, channel_index_start,
                                   channel_index_end,
                                   v_residual[time_index].data());
  }

  // Subtract all directions with their current solutions
  for (size_t d = 0; d != NDirections(); ++d)
    AddOrSubtractDirection<false>(solver_buffer, v_residual,
                                  channel_block_index, d, solutions);

  const std::vector<std::vector<Complex>> v_copy = v_residual;

  for (size_t direction = 0; direction != NDirections(); ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    AddOrSubtractDirection<true>(solver_buffer, v_residual, channel_block_index,
                                 direction, solutions);

    SolveDirection(solver_buffer, channel_block_index, v_residual, direction,
                   solutions, next_solutions);
  }
}

void IterativeDiagonalSolver::SolveDirection(
    const SolverBuffer& solver_buffer, size_t channel_block_index,
    const std::vector<std::vector<Complex>>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  const size_t channel_index_start =
      channel_block_index * n_channels_ / n_channel_blocks_;
  const size_t channel_index_end =
      (channel_block_index + 1) * n_channels_ / n_channel_blocks_;
  const size_t n_times = solver_buffer.NTimes();

  // Now calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^H
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  std::vector<aocommon::MC2x2Diag> numerator(NAntennas(),
                                             aocommon::MC2x2Diag::Zero());
  std::vector<double> denominator(NAntennas() * 2, 0.0);

  // Iterate over all data
  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    const Complex* data_ptr = v_residual[time_index].data();
    for (size_t baseline = 0; baseline != ant1_.size(); ++baseline) {
      const size_t antenna1 = ant1_[baseline];
      const size_t antenna2 = ant2_[baseline];
      if (antenna1 != antenna2) {
        const DComplex* solution_ant1 =
            &solutions[(antenna1 * NDirections() + direction) * 2];
        const DComplex* solution_ant2 =
            &solutions[(antenna2 * NDirections() + direction) * 2];
        const Complex* model_ptr = solver_buffer.ModelDataPointer(
            time_index, direction, baseline, channel_index_start);
        for (size_t ch = channel_index_start; ch != channel_index_end; ++ch) {
          // Calculate the contribution of this baseline for antenna1
          const aocommon::MC2x2 data(data_ptr);
          const aocommon::MC2x2Diag solution1(solution_ant2[0],
                                              solution_ant2[1]);
          const aocommon::MC2x2 cor_model_transp1(
              solution1 * HermTranspose(aocommon::MC2x2(model_ptr)));
          numerator[antenna1] += Diagonal(data * cor_model_transp1);
          // The indices (0, 2 / 1, 3) are following from the fact that we want
          // the contribution of antenna2's "X" polarization, and the matrix is
          // ordered [ XX XY / YX YY ].
          denominator[antenna1 * 2] +=
              std::norm(cor_model_transp1[0]) + std::norm(cor_model_transp1[2]);
          denominator[antenna1 * 2 + 1] +=
              std::norm(cor_model_transp1[1]) + std::norm(cor_model_transp1[3]);

          // Calculate the contribution of this baseline for antenna2
          // data_ba = data_ab^H, etc., therefore, numerator and denominator
          // become:
          // - num = data_ab^H * solutions_a * model_ab
          // - den = norm(model_ab^H * solutions_a)
          const aocommon::MC2x2Diag solution2(solution_ant1[0],
                                              solution_ant1[1]);
          const aocommon::MC2x2 cor_model2(solution2 *
                                           aocommon::MC2x2(model_ptr));

          numerator[antenna2] += Diagonal(HermTranspose(data) * cor_model2);
          denominator[antenna2 * 2] +=
              std::norm(cor_model2[0]) + std::norm(cor_model2[2]);
          denominator[antenna2 * 2 + 1] +=
              std::norm(cor_model2[1]) + std::norm(cor_model2[3]);

          data_ptr += 4;  // Skip to next 2x2 matrix
          model_ptr += 4;
        }
      } else {
        // skip autocorrelation, therefore skip nr. channels data matrices:
        data_ptr += 4 * (channel_index_end - channel_index_start);
      }
    }
  }
  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    DComplex* destination =
        &next_solutions[(ant * NDirections() + direction) * 2];

    if (denominator[ant * 2] == 0.0)
      destination[0] = std::numeric_limits<double>::quiet_NaN();
    else
      destination[0] = numerator[ant][0] / denominator[ant * 2];

    if (denominator[ant * 2 + 1] == 0.0)
      destination[1] = std::numeric_limits<double>::quiet_NaN();
    else
      destination[1] = numerator[ant][1] / denominator[ant * 2 + 1];
  }
}

template <bool Add>
void IterativeDiagonalSolver::AddOrSubtractDirection(
    const SolverBuffer& solver_buffer,
    std::vector<std::vector<Complex>>& v_residual, size_t channel_block_index,
    size_t direction, const std::vector<DComplex>& solutions) {
  const size_t channel_index_start =
      channel_block_index * n_channels_ / n_channel_blocks_;
  const size_t channel_index_end =
      (channel_block_index + 1) * n_channels_ / n_channel_blocks_;
  const size_t n_times = solver_buffer.NTimes();

  for (size_t time_index = 0; time_index != n_times; ++time_index) {
    Complex* data_ptr = v_residual[time_index].data();
    for (size_t baseline = 0; baseline != ant1_.size(); ++baseline) {
      const size_t antenna1 = ant1_[baseline];
      const size_t antenna2 = ant2_[baseline];
      const Complex* model_ptr = solver_buffer.ModelDataPointer(
          time_index, direction, baseline, channel_index_start);
      for (size_t ch = channel_index_start; ch != channel_index_end; ++ch) {
        const DComplex* solution1 =
            &solutions[(antenna1 * NDirections() + direction) * 2];
        const DComplex* solution2 =
            &solutions[(antenna2 * NDirections() + direction) * 2];

        if (Add) {
          data_ptr[0] += Complex(solution1[0]) * model_ptr[0] *
                         std::conj(Complex(solution2[0]));
          data_ptr[1] += Complex(solution1[0]) * model_ptr[1] *
                         std::conj(Complex(solution2[1]));
          data_ptr[2] += Complex(solution1[1]) * model_ptr[2] *
                         std::conj(Complex(solution2[0]));
          data_ptr[3] += Complex(solution1[1]) * model_ptr[3] *
                         std::conj(Complex(solution2[1]));
        } else {
          data_ptr[0] -= Complex(solution1[0]) * model_ptr[0] *
                         std::conj(Complex(solution2[0]));
          data_ptr[1] -= Complex(solution1[0]) * model_ptr[1] *
                         std::conj(Complex(solution2[1]));
          data_ptr[2] -= Complex(solution1[1]) * model_ptr[2] *
                         std::conj(Complex(solution2[0]));
          data_ptr[3] -= Complex(solution1[1]) * model_ptr[3] *
                         std::conj(Complex(solution2[1]));
        }

        model_ptr += 4;  // Goto the next 2x2 matrix
        data_ptr += 4;
      }
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
