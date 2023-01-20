// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeScalarSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/parallelfor.h>

#include <algorithm>

using aocommon::MC2x2;
using aocommon::MC2x2F;

namespace dp3 {
namespace ddecal {

IterativeScalarSolver::SolveResult IterativeScalarSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(NChannelBlocks());

  SolveResult result;

  // Visibility vector v_residual[cb][vis] of size NChannelBlocks() x
  // n_visibilities
  std::vector<std::vector<MC2x2F>> v_residual(NChannelBlocks());
  // The following loop allocates all structures
  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    next_solutions[ch_block].resize(NSolutions() * NAntennas());
    v_residual[ch_block].resize(data.ChannelBlock(ch_block).NVisibilities());
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
    MakeSolutionsFinite1Pol(solutions);

    aocommon::ParallelFor<size_t> loop(GetNThreads());
    loop.Run(0, NChannelBlocks(),
             [&](size_t ch_block, [[maybe_unused]] size_t thread) {
               PerformIteration(data.ChannelBlock(ch_block),
                                v_residual[ch_block], solutions[ch_block],
                                next_solutions[ch_block]);
             });

    Step(solutions, next_solutions);

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes);
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

void IterativeScalarSolver::PerformIteration(
    const SolveData::ChannelBlockData& cb_data, std::vector<MC2x2F>& v_residual,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction)
    AddOrSubtractDirection<false>(cb_data, v_residual, direction, solutions);

  const std::vector<MC2x2F> v_copy = v_residual;

  for (size_t direction = 0; direction != NDirections(); ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    AddOrSubtractDirection<true>(cb_data, v_residual, direction, solutions);

    SolveDirection(cb_data, v_residual, direction, solutions, next_solutions);
  }
}

void IterativeScalarSolver::SolveDirection(
    const SolveData::ChannelBlockData& cb_data,
    const std::vector<MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  const size_t n_dir_solutions = cb_data.NSolutionsForDirection(direction);
  std::vector<std::complex<double>> numerator(NAntennas() * n_dir_solutions,
                                              0.0);
  std::vector<double> denominator(NAntennas() * n_dir_solutions, 0.0);

  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  const std::vector<uint32_t>& solution_map = cb_data.SolutionMap(direction);
  const std::vector<MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const size_t antenna_1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna_2 = cb_data.Antenna2Index(vis_index);
    const uint32_t solution_index = solution_map[vis_index];
    const Complex solution_ant_1(
        solutions[antenna_1 * NSolutions() + solution_index]);
    const Complex solution_ant_2(
        solutions[antenna_2 * NSolutions() + solution_index]);
    const MC2x2F& data = v_residual[vis_index];
    const MC2x2F& model = model_vector[vis_index];

    const uint32_t rel_solution_index = solution_index - solution_map[0];
    // Calculate the contribution of this baseline for antenna_1
    const MC2x2F cor_model_herm_1(HermTranspose(model) * solution_ant_2);
    const size_t full_solution_1_index =
        antenna_1 * n_dir_solutions + rel_solution_index;
    numerator[full_solution_1_index] += Trace(data * cor_model_herm_1);
    denominator[full_solution_1_index] += Norm(cor_model_herm_1);

    // Calculate the contribution of this baseline for antenna2
    const MC2x2F cor_model_2(model * solution_ant_1);
    const size_t full_solution_2_index =
        antenna_2 * n_dir_solutions + rel_solution_index;
    numerator[full_solution_2_index] +=
        Trace(HermTranspose(data) * cor_model_2);
    denominator[full_solution_2_index] += Norm(cor_model_2);
  }

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    for (size_t rel_sol = 0; rel_sol != n_dir_solutions; ++rel_sol) {
      const uint32_t solution_index = rel_sol + solution_map[0];
      DComplex& destination =
          next_solutions[ant * NSolutions() + solution_index];
      const size_t index = ant * n_dir_solutions + rel_sol;
      if (denominator[index] == 0.0)
        destination = std::numeric_limits<float>::quiet_NaN();
      else
        destination = numerator[index] / denominator[index];
    }
  }
}

template <bool Add>
void IterativeScalarSolver::AddOrSubtractDirection(
    const SolveData::ChannelBlockData& cb_data, std::vector<MC2x2F>& v_residual,
    size_t direction, const std::vector<DComplex>& solutions) {
  const std::vector<MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  const size_t n_visibilities = cb_data.NVisibilities();
  const std::vector<uint32_t>& solution_map = cb_data.SolutionMap(direction);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
    const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
    const uint32_t solution_index = solution_map[vis_index];
    const Complex solution_1(
        solutions[antenna_1 * NSolutions() + solution_index]);
    const Complex solution_2_conj = std::conj(
        Complex(solutions[antenna_2 * NSolutions() + solution_index]));
    MC2x2F& data = v_residual[vis_index];
    const MC2x2F& model = model_vector[vis_index];
    if (Add) {
      data.AddWithFactorAndAssign(model, solution_1 * solution_2_conj);
    } else {
      data.AddWithFactorAndAssign(model, -solution_1 * solution_2_conj);
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
